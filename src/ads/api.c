#include "ads/api.h"
#include "ads/isax_index.h"
#include "ads/inmemory_index_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/sax/ts.h"
#include "ads/sfa/sfa.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/spartan.h"
#include <stdlib.h>
#include <float.h>
#include <stdint.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>

struct messi_index {
    isax_index *index;
};

static void populate_root_nodes(isax_index *index, node_list *list);
static void prepare_sfa_bins_if_needed(isax_index *index, const char *path, long ts_num);
static void finalize_sfa_bins_if_needed(isax_index *index);
static void prepare_spartan_bins_if_needed(isax_index *index, const char *path, long ts_num);
static void finalize_spartan_bins_if_needed(isax_index *index);

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
        params->n_coefficients);

    if (settings == NULL) {
        free(wrapper);
        return NULL;
    }

    wrapper->index = isax_index_init(settings);
    if (wrapper->index == NULL) {
        free(wrapper);
        return NULL;
    }

    return wrapper;
}

void messi_index_destroy(messi_index *index) {
    if (index == NULL) {
        return;
    }
    if (index->index != NULL) {
        if (index->index->settings && index->index->settings->function_type == 4) {
            sfa_free_bins(index->index);
        }
        if (index->index->settings && index->index->settings->function_type == 5) {
            spartan_free_bins(index->index);
        }
        MESSI2_index_destroy(index->index, index->index->first_node);
    }
    free(index);
}

int messi_index_add_file(messi_index *index, const char *path, long ts_num) {
    if (index == NULL || index->index == NULL || path == NULL) {
        return -1;
    }
    prepare_sfa_bins_if_needed(index->index, path, ts_num);
    prepare_spartan_bins_if_needed(index->index, path, ts_num);
    index_creation_pRecBuf(path, ts_num, 0, 0, index->index);
    finalize_sfa_bins_if_needed(index->index);
    finalize_spartan_bins_if_needed(index->index);

    return 0;
}

int messi_index_search(messi_index *index,
                       const float *queries,
                       size_t nq,
                       size_t dim,
                       size_t k,
                       float *distances,
                       long *labels) {
    if (index == NULL || index->index == NULL || queries == NULL || distances == NULL || labels == NULL) {
        return -1;
    }
    if (dim != (size_t) index->index->settings->timeseries_size) {
        return -2;
    }

    ts_type *paa_buffer = malloc(sizeof(ts_type) * index->index->settings->n_segments);

    fftw_workspace fftw = {0};

    unsigned long ts_length = index->index->settings->timeseries_size;
    if (index->index->settings->function_type == 4) {
        fftw_workspace_init(&fftw, ts_length);
    }

    for (size_t i = 0; i < nq; ++i) {
        node_list nlist = {.nlist = NULL, .node_amount = 0};
        populate_root_nodes(index->index, &nlist);
        ts_type *ts = (ts_type *) (queries + i * dim);

        if (index->index->settings->function_type == 4) {
            //SFA: parse ts and make fft representation
            memcpy(fftw.ts, ts, sizeof(ts_type) * ts_length);

            int use_best = index->index->settings->n_coefficients != 0;
            fft_from_ts(
				index->index,
				index->index->settings->n_segments,
				use_best, &fftw);

            memcpy(paa_buffer, fftw.transform, sizeof(ts_type) * index->index->settings->n_segments);
        } else if (index->index->settings->function_type == 5) {
            pca_from_ts(index->index, ts, paa_buffer);
        } else {
            paa_from_ts(ts,
                        paa_buffer,
                        index->index->settings);
        }
        query_result res = exact_search_MESSI(
			(ts_type *) (queries + i * dim),
            paa_buffer,
            index->index,
            &nlist,
            FLT_MAX,
            -1);
        distances[i] = res.distance;
        labels[i] = res.node ? (long) (intptr_t) res.node : -1;
        if (nlist.nlist != NULL) {
            free(nlist.nlist);
        }
    }
    if (index->index->settings->function_type == 4) {
        fftw_workspace_destroy(&fftw);
    }
    free(paa_buffer);
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

static void prepare_sfa_bins_if_needed(isax_index *index, const char *path, long ts_num) {
    if (index == NULL || index->settings->function_type != 4) {
        return;
    }
    if (sfa_bins_init(index) != SUCCESS) {
        fprintf(stderr, "warning: failed to initialize SFA bins.\n");
        return;
    }
    sfa_set_bins(index, path, ts_num, maxquerythread, 1, !index->settings->is_norm);
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

static void prepare_spartan_bins_if_needed(isax_index *index, const char *path, long ts_num) {
    if (index == NULL || index->settings->function_type != 5) {
        return;
    }
    if (spartan_bins_init(index) != SUCCESS) {
        fprintf(stderr, "warning: failed to initialize SPARTAN bins.\n");
        return;
    }
    spartan_set_bins(index, path, ts_num, maxquerythread, 1, !index->settings->is_norm);
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
