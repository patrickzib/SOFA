#include "ads/api.h"
#include "ads/isax_index.h"
#include "ads/inmemory_index_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/sax/ts.h"
#include "ads/calc_utils.h"
#include "ads/sfa/sfa.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"
#include "ads/sffa/sffa.h"
#include "ads/trie/trie.h"
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* The CLI owns this global in executable builds.  The standalone API links
 * query engines directly, so it provides the same default independently. */
int query_report_interval = 10;

struct messi_index {
    isax_index *index;
    int filetype_int;
    int built;
};

static int fanout_to_bits(int fanout) {
    int bits = 0;
    if (fanout < 2 || fanout > 256 || (fanout & (fanout - 1)) != 0) return -1;
    while (fanout > 1) {
        ++bits;
        fanout >>= 1;
    }
    return bits;
}

static void populate_root_nodes(isax_index *index, node_list *list);
static enum response prepare_sfa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int);
static void finalize_sfa_bins_if_needed(isax_index *index);
static enum response prepare_spartan_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int);
static void finalize_spartan_bins_if_needed(isax_index *index);
static enum response prepare_pisa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int);
static void finalize_pisa_bins_if_needed(isax_index *index);
static enum response prepare_sffa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int);
static void finalize_sffa_bins_if_needed(isax_index *index);

static messi_index *messi_alloc(void) {
    messi_index *idx = (messi_index *) calloc(1, sizeof(messi_index));
    return idx;
}

messi_index *messi_index_create(const messi_index_params *params) {
    if (params == NULL) {
        return NULL;
    }

    int requested_threads = params->max_query_threads;
    if (requested_threads <= 0) {
        requested_threads = 1;
    }
    maxquerythread = requested_threads;

    int requested_queues = params->queue_count;
    if (requested_queues <= 0) {
        requested_queues = requested_threads;
    }
    N_PQUEUE = requested_queues;

    messi_index *wrapper = messi_alloc();
    if (wrapper == NULL) {
        return NULL;
    }

    char cwd_buffer[PATH_MAX];
    const char *root_directory = params->root_directory;
    if (root_directory == NULL || *root_directory == '\0') {
        if (getcwd(cwd_buffer, sizeof(cwd_buffer)) != NULL) {
            root_directory = cwd_buffer;
        } else {
            root_directory = ".";
        }
    }

    const messi_index_type index_type = params->index_type == MESSI_INDEX_TRIE
                                            ? MESSI_INDEX_TRIE : MESSI_INDEX_ISAX;
    isax_index_settings *settings = isax_index_settings_init(
        root_directory,
        params->timeseries_size,
        params->n_segments,
        params->sax_bit_cardinality,
        params->max_leaf_size,
        params->min_leaf_size,
        params->initial_leaf_buffer_size,
        params->max_total_buffer_size,
        params->initial_fbl_buffer_size,
        params->total_loaded_leaves,
        params->tight_bound,
        params->aggressive_check,
        1,
        params->function_type,
        1,
        params->simd,
        params->sample_size,
        params->is_norm,
        params->histogram_type,
        params->sample_type,
        params->n_coefficients,
        index_type);

    if (settings == NULL) {
        free(wrapper);
        return NULL;
    }

    settings->sampling_seed = params->sampling_seed == 0 ? 1 : params->sampling_seed;
    settings->sffa_order = params->sffa_order;
    settings->sffa_auto_order = params->sffa_auto_order;
    if (!isfinite(settings->sffa_order)) {
        free(settings);
        free(wrapper);
        return NULL;
    }
    settings->node_split_criterion = params->node_split_criterion == 0
                                        ? 1 : params->node_split_criterion;
    /* Match the CLI/default constructor: node MBRs remain enabled for iSAX. */
    settings->isax_node_mbr = 1;
    settings->isax_record_mbr_suffix_bound = params->isax_record_mbr_suffix_bound;
    settings->isax_record_lb_table = params->isax_record_lb_table;
    settings->isax_mbr_dimensions = params->isax_mbr_dimensions == 0
                                        ? 32 : params->isax_mbr_dimensions;
    settings->trie_bound_dimensions = params->trie_bound_dimensions;
    settings->trie_split_dimensions = params->trie_split_dimensions;
    settings->trie_record_mbr_suffix_bound = params->trie_record_mbr_suffix_bound;
    settings->trie_leaf_ivf = params->trie_leaf_ivf;
    settings->trie_fanout = params->trie_fanout == 0 ? 8 : params->trie_fanout;
    settings->trie_dynamic_alphabet = params->trie_dynamic_alphabet;
    settings->trie_min_bits = fanout_to_bits(params->trie_min_fanout == 0
                                                  ? 2 : params->trie_min_fanout);
    settings->trie_max_bits = fanout_to_bits(params->trie_max_fanout == 0
                                                  ? 16 : params->trie_max_fanout);
    settings->trie_alphabet_budget_bits = params->trie_alphabet_budget_bits == 0
                                                ? 3 : params->trie_alphabet_budget_bits;
    settings->dynamic_root_split_variance =
        index_type == MESSI_INDEX_ISAX && params->dynamic_root_split_variance;
    if (settings->node_split_criterion < 1 || settings->node_split_criterion > 4 ||
        settings->trie_min_bits < 1 || settings->trie_max_bits < settings->trie_min_bits ||
        settings->trie_max_bits > 8 ||
        settings->trie_alphabet_budget_bits < settings->trie_min_bits ||
        settings->trie_alphabet_budget_bits > settings->trie_max_bits) {
        free(settings);
        free(wrapper);
        return NULL;
    }

    wrapper->index = isax_index_init_inmemory(settings);
    if (wrapper->index == NULL) {
        free(wrapper);
        return NULL;
    }
    int filetype_int = params->filetype_int;
    if (filetype_int != 0 && filetype_int != 1) {
        filetype_int = 0;
    }
    wrapper->filetype_int = filetype_int;

    return wrapper;
}

void messi_index_destroy(messi_index *index) {
    if (index == NULL) {
        return;
    }
    if (index->index != NULL) {
        if (index->index->settings && index->index->bins != NULL &&
            index->index->settings->function_type == 4) {
            sfa_free_bins(index->index);
        }
        if (index->index->settings && index->index->settings->function_type == 5) {
            spartan_free_bins(index->index);
        }
        if (index->index->settings && index->index->bins != NULL &&
            index->index->settings->function_type == 6) {
            pisa_free_bins(index->index);
        }
        if (index->index->settings && index->index->bins != NULL &&
            index->index->settings->function_type == 7) {
            sffa_free_bins(index->index);
        }
        if (index->built && index->index->settings &&
            index->index->settings->index_type == MESSI_INDEX_TRIE) {
            symbolic_trie_destroy(index->index);
        }
        /* rawfile is the process-global in-memory raw-series buffer used by
         * both pRecBuf and symbolic-trie exact refinement. */
        if (index->built && rawfile != NULL) {
            free(rawfile);
            rawfile = NULL;
        }
        isax_index_pRecBuf_destroy(index->index, NULL, maxquerythread);
    }
    free(index);
}

int messi_index_add_file(messi_index *index, const char *path, long ts_num, int dynamic_index) {
    if (index == NULL || index->index == NULL || path == NULL || ts_num <= 0 || index->built) {
        return -1;
    }
    int apply_znorm = index->index->settings ? !index->index->settings->is_norm : 0;
    if (prepare_sfa_bins_if_needed(index->index, path, ts_num, index->filetype_int) != SUCCESS ||
        prepare_spartan_bins_if_needed(index->index, path, ts_num, index->filetype_int) != SUCCESS ||
        prepare_pisa_bins_if_needed(index->index, path, ts_num, index->filetype_int) != SUCCESS ||
        prepare_sffa_bins_if_needed(index->index, path, ts_num, index->filetype_int) != SUCCESS) {
        return -3;
    }
    enum response build_result;
    if (index->index->settings->index_type == MESSI_INDEX_TRIE) {
        build_result = symbolic_trie_build(index->index, path, ts_num,
                                           index->filetype_int, apply_znorm);
    } else {
        index_creation_pRecBuf(path, ts_num, index->filetype_int, apply_znorm,
                               index->index, dynamic_index);
        build_result = SUCCESS;
    }
    if (build_result != SUCCESS) {
        finalize_sfa_bins_if_needed(index->index);
        finalize_spartan_bins_if_needed(index->index);
        finalize_pisa_bins_if_needed(index->index);
        finalize_sffa_bins_if_needed(index->index);
        return -4;
    }
    finalize_sfa_bins_if_needed(index->index);
    finalize_spartan_bins_if_needed(index->index);
    finalize_pisa_bins_if_needed(index->index);
    finalize_sffa_bins_if_needed(index->index);
    index->built = 1;
    return 0;
}

double messi_index_sffa_order(const messi_index *index) {
    if (index == NULL || index->index == NULL || index->index->settings == NULL ||
        index->index->settings->function_type != MESSI_TRANSFORM_SFFA) {
        return NAN;
    }
    return index->index->settings->sffa_order;
}

int messi_index_search(messi_index *index,
                       const float *queries,
                       size_t nq,
                       size_t dim,
                       size_t k,
                       float *distances,
                       long *labels,
                       int dynamic_index) {
    if (index == NULL || index->index == NULL || queries == NULL || distances == NULL || labels == NULL) {
        return -1;
    }
    if (dim != (size_t) index->index->settings->timeseries_size) {
        return -2;
    }

    ts_type *paa_buffer = malloc(sizeof(ts_type) * index->index->settings->n_segments);
    ts_type *normalized_query = NULL;
    if (paa_buffer == NULL) return -3;
    if (!index->index->settings->is_norm) {
        normalized_query = malloc(sizeof(*normalized_query) * dim);
        if (normalized_query == NULL) { free(paa_buffer); return -3; }
    }

    fftw_workspace fftw = {0};
    sffa_workspace sffa = {0};

    unsigned long ts_length = index->index->settings->timeseries_size;
    if (index->index->settings->function_type == 4 || index->index->settings->function_type == 6) {
        fftw_workspace_init(&fftw, ts_length);
    }
    if (index->index->settings->function_type == 7) {
        sffa_workspace_init(&sffa, ts_length);
    }

    for (size_t i = 0; i < nq; ++i) {
        node_list nlist = {.nlist = NULL, .node_amount = 0};
        populate_root_nodes(index->index, &nlist);
        ts_type *ts = (ts_type *) (queries + i * dim);
        if (normalized_query != NULL) {
            memcpy(normalized_query, ts, sizeof(*normalized_query) * dim);
            znorm(normalized_query, (int) dim);
            ts = normalized_query;
        }

        if (index->index->settings->function_type == 4) {
            //SFA: parse ts and make fft representation
            memcpy(fftw.ts, ts, sizeof(ts_type) * ts_length);

            int use_best = index->index->settings->n_coefficients != 0;
            if (fft_from_ts(index->index, index->index->settings->n_segments, use_best,
                            &fftw) != SUCCESS) {
                free(nlist.nlist);
                fftw_workspace_destroy(&fftw);
                free(normalized_query);
                free(paa_buffer);
                return -5;
            }

            memcpy(paa_buffer, fftw.transform, sizeof(ts_type) * index->index->settings->n_segments);
        } else if (index->index->settings->function_type == 5) {
            pca_from_ts(index->index, ts, paa_buffer);
        } else if (index->index->settings->function_type == 6) {
            if (pisa_pca_from_ts(index->index, ts, paa_buffer, &fftw) != SUCCESS) {
                free(nlist.nlist);
                fftw_workspace_destroy(&fftw);
                free(normalized_query);
                free(paa_buffer);
                return -5;
            }
        } else if (index->index->settings->function_type == 7) {
            if (sffa_project(index->index, ts, paa_buffer, &sffa) != SUCCESS) {
                free(nlist.nlist);
                sffa_workspace_destroy(&sffa);
                free(normalized_query);
                free(paa_buffer);
                return -5;
            }
        } else {
            paa_from_ts(ts,
                        paa_buffer,
                        index->index->settings);
        }
        query_result res;
        if (index->index->settings->index_type == MESSI_INDEX_TRIE) {
            res = symbolic_trie_exact_search(index->index, ts, paa_buffer, FLT_MAX);
        } else {
            res = exact_search_MESSI(
                ts, paa_buffer, index->index,
                &nlist, FLT_MAX, -1, dynamic_index);
        }
        distances[i] = res.distance;
        labels[i] = res.record_position == QUERY_RESULT_NO_POSITION
                        ? -1
                        : (long) (res.record_position / index->index->settings->timeseries_size);
        if (nlist.nlist != NULL) {
            free(nlist.nlist);
        }
    }
    if (index->index->settings->function_type == 4 || index->index->settings->function_type == 6) {
        fftw_workspace_destroy(&fftw);
    }
    if (index->index->settings->function_type == 7) sffa_workspace_destroy(&sffa);
    free(normalized_query);
    free(paa_buffer);
    return 0;
}

int messi_index_pca_transform(messi_index *index,
                              const float *queries,
                              size_t nq,
                              size_t dim,
                              float *out,
                              size_t out_dim) {
    if (index == NULL || index->index == NULL || queries == NULL || out == NULL) {
        return -1;
    }
    if (index->index->settings == NULL) {
        return -1;
    }
    if (dim != (size_t) index->index->settings->timeseries_size) {
        return -2;
    }
    if (out_dim != (size_t) index->index->settings->n_segments) {
        return -3;
    }
    if (index->index->settings->function_type != 5) {
        return -4;
    }

    for (size_t i = 0; i < nq; ++i) {
        const ts_type *ts = (const ts_type *) (queries + i * dim);
        ts_type *row_out = (ts_type *) (out + i * out_dim);
        if (pca_from_ts(index->index, ts, row_out) != SUCCESS) {
            return -5;
        }
    }
    return 0;
}
static void populate_root_nodes(isax_index *index, node_list *list) {
    if (list == NULL || index == NULL) {
        return;
    }
    size_t capacity = index->settings->root_nodes_size;
    list->nlist = malloc(sizeof(isax_node *) * capacity);
    list->node_amount = 0;
    if (list->nlist == NULL) {
        return;
    }
    isax_node *current_root_node = index->first_node;
    while (current_root_node != NULL && (size_t) list->node_amount < capacity) {
        list->nlist[list->node_amount++] = current_root_node;
        current_root_node = current_root_node->next;
    }
}

static enum response prepare_sfa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int) {
    if (index == NULL || index->settings->function_type != 4) {
        return SUCCESS;
    }
    if (sfa_bins_init(index) != SUCCESS) {
        fprintf(stderr, "error: failed to initialize SFA bins.\n");
        return FAILURE;
    }
    return sfa_set_bins(index, path, ts_num, maxquerythread, filetype_int, !index->settings->is_norm);
}

static void finalize_sfa_bins_if_needed(isax_index *index) {
    if (index == NULL || index->settings->function_type != 4 || index->bins == NULL || index->binsv == NULL) {
        return;
    }
    int n_segments = index->settings->n_segments;
    int slice = index->settings->sax_alphabet_cardinality - 1;
    if (slice <= 0) {
        return;
    }
    for (int i = 0; i < n_segments; i++) {
        memcpy(&index->binsv[i * slice], index->bins[i], sizeof(ts_type) * slice);
    }
}

static enum response prepare_spartan_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int) {
    if (index == NULL || index->settings->function_type != 5) {
        return SUCCESS;
    }
    if (spartan_bins_init(index) != SUCCESS) {
        fprintf(stderr, "error: failed to initialize SPARTAN bins.\n");
        return FAILURE;
    }
    return spartan_set_bins(index, path, ts_num, maxquerythread, filetype_int, !index->settings->is_norm);
}

static void finalize_spartan_bins_if_needed(isax_index *index) {
    if (index == NULL || index->settings->function_type != 5 || index->bins == NULL || index->binsv == NULL) {
        return;
    }
    int n_segments = index->settings->n_segments;
    int slice = index->settings->sax_alphabet_cardinality - 1;
    if (slice <= 0) {
        return;
    }
    for (int i = 0; i < n_segments; i++) {
        memcpy(&index->binsv[i * slice], index->bins[i], sizeof(ts_type) * slice);
    }
}

static enum response prepare_pisa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int) {
    if (index == NULL || index->settings->function_type != 6) {
        return SUCCESS;
    }
    if (pisa_bins_init(index) != SUCCESS) {
        fprintf(stderr, "error: failed to initialize PISA bins.\n");
        return FAILURE;
    }
    return pisa_set_bins(index, path, ts_num, maxquerythread, filetype_int, !index->settings->is_norm);
}

static void finalize_pisa_bins_if_needed(isax_index *index) {
    if (index == NULL || index->settings->function_type != 6 || index->bins == NULL || index->binsv == NULL) {
        return;
    }
    int n_segments = index->settings->n_segments;
    int slice = index->settings->sax_alphabet_cardinality - 1;
    if (slice <= 0) {
        return;
    }
    for (int i = 0; i < n_segments; i++) {
        memcpy(&index->binsv[i * slice], index->bins[i], sizeof(ts_type) * slice);
    }
}

static enum response prepare_sffa_bins_if_needed(isax_index *index, const char *path, long ts_num, int filetype_int) {
    if (index == NULL || index->settings->function_type != 7) return SUCCESS;
    if (sffa_bins_init(index) != SUCCESS) {
        fprintf(stderr, "error: failed to initialize SFFA bins.\n");
        return FAILURE;
    }
    return sffa_set_bins(index, path, ts_num, maxquerythread, filetype_int, !index->settings->is_norm);
}

static void finalize_sffa_bins_if_needed(isax_index *index) {
    (void) index;
    /* SFFA's bin rows alias its contiguous binsv allocation. */
}
