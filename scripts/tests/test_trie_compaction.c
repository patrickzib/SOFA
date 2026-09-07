/* Include the implementation to test private arena ownership and layout without
 * adding a public inspection API. Link this test with libads and its usual libs. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>

int query_report_interval = 0;

static int fail_next_allocation;
static void *compaction_test_malloc(size_t size) {
    if (fail_next_allocation) { fail_next_allocation = 0; return NULL; }
    return malloc(size);
}
#define malloc compaction_test_malloc
#include "../../src/ads/trie/trie.c"
#undef malloc

static int check_compaction(int dimensions, int split) {
    const int order[] = {4, 0, 6, 1, 5, 2, 3};
    struct symbolic_trie_index trie = {0};
    trie.dimensions = dimensions;
    trie.word_arena = malloc(7 * (size_t) dimensions);
    trie.root = trie_node_create(dimensions, NULL, (trie_dimension_mask) {0, 0});
    if (!trie.word_arena || !trie.root) return 0;
    for (int r = 0; r < 7; ++r)
        for (int d = 0; d < dimensions; ++d)
            trie.word_arena[r * dimensions + d] = (sax_type) (r * 19 + d);
    symbolic_trie_node *leaves[3] = {trie.root, NULL, NULL};
    int count = 1;
    if (split) {
        trie.root->leaf = 0;
        trie.root->split_fanout = 4;
        count = 3;
        for (int i = 0; i < count; ++i) {
            leaves[i] = trie_node_create(dimensions, trie.root, (trie_dimension_mask) {0, 0});
            if (!leaves[i]) return 0;
            trie.root->children[i + 1] = leaves[i]; /* Also cover a missing child. */
        }
    }
    int record = 0;
    for (int i = 0; i < count; ++i) {
        const int size = split ? (i == 0 ? 3 : i == 1 ? 0 : 4) : 7;
        for (int j = 0; j < size; ++j, ++record)
            if (!trie_leaf_append(leaves[i], trie.word_arena + order[record] * dimensions,
                                  order[record] * 64, dimensions)) return 0;
        if (!size) continue;
        leaves[i]->record_raw_radii = malloc(size * sizeof(float));
        leaves[i]->clusters = malloc(sizeof(trie_leaf_cluster));
        if (!leaves[i]->record_raw_radii || !leaves[i]->clusters) return 0;
        leaves[i]->clusters[0] = (trie_leaf_cluster) {0, size};
        leaves[i]->cluster_count = 1;
        for (int j = 0; j < size; ++j)
            leaves[i]->record_raw_radii[j] = leaves[i]->positions[j] / 64 + 0.5f;
    }
    sax_type *original = trie.word_arena;
    sax_type *first_word = leaves[0]->words[0];
    fail_next_allocation = 1;
    int ok = !trie_compact_words(&trie, 7) && trie.word_arena == original &&
             leaves[0]->words[0] == first_word;
    if (dimensions > 1 && trie_compact_words(&trie, SIZE_MAX)) ok = 0;
    for (int pass = 0; pass < 2 && ok; ++pass) {
        if (!trie_compact_words(&trie, 7)) { ok = 0; break; }
        record = 0;
        for (int i = 0; i < count; ++i) {
            symbolic_trie_node *leaf = leaves[i];
            if (leaf->size && (leaf->cluster_count != 1 || leaf->clusters[0].offset != 0 ||
                              leaf->clusters[0].size != leaf->size)) ok = 0;
            for (int j = 0; j < leaf->size; ++j, ++record) {
                if (leaf->words[j] != trie.word_arena + record * dimensions ||
                    leaf->positions[j] != (file_position_type) order[record] * 64 ||
                    leaf->record_raw_radii[j] != order[record] + 0.5f) ok = 0;
                for (int d = 0; d < dimensions; ++d)
                    if (leaf->words[j][d] != (sax_type) (order[record] * 19 + d) ||
                        leaf->words[j][d] < leaf->min_word[d] ||
                        leaf->words[j][d] > leaf->max_word[d]) ok = 0;
            }
        }
        if (record != 7) ok = 0;
    }
    trie_node_destroy(trie.root);
    free(trie.word_arena);
    return ok;
}

/* A winning seed is skipped during traversal, but its record ID must survive
 * both the serial and parallel API paths, including an equal supplied BSF. */
static int check_seed_search(void) {
    enum { D = 16 };
    isax_index index = {0};
    isax_index_settings settings = {0};
    struct symbolic_trie_index trie = {0};
    float boundaries[D][255], *bins[D], query[D];
    sax_type cards[D], words[2][D], query_word[D];
    settings.n_segments = settings.timeseries_size = D;
    settings.sax_bit_cardinality = 8;
    settings.sax_alphabet_cardinality = 256;
    settings.function_type = 5;
    settings.SIMD_flag = 1;
    settings.trie_streaming_leaf_scan = 1;
    settings.max_sax_cardinalities = cards;
    index.settings = &settings;
    index.bins = bins;
    index.trie = &trie;
    trie.dimensions = trie.bound_dimensions = D;
    trie.root = trie_node_create(D, NULL, (trie_dimension_mask) {0, 0});
    if (!trie.root) return 0;
    trie.root->leaf = 0;
    trie.root->split_dimension = 0;
    trie.root->split_fanout = 2;
    rawfile = malloc(2 * D * sizeof(*rawfile));
    if (!rawfile) return 0;
    for (int d = 0; d < D; ++d) {
        bins[d] = boundaries[d]; cards[d] = 8; query[d] = 1.0f;
        words[0][d] = 64; words[1][d] = 144;
        rawfile[d] = -4.0f; rawfile[D + d] = 1.0f;
        for (int s = 0; s < 255; ++s) boundaries[d][s] = (s - 127.0f) / 16.0f;
    }
    for (int r = 0; r < 2; ++r) {
        trie.root->children[r] = trie_node_create(D, trie.root, (trie_dimension_mask) {0, 0});
        if (!trie.root->children[r] ||
            !trie_leaf_append(trie.root->children[r], words[r], r * D, D)) return 0;
        trie_node_update_mbb(trie.root, words[r], D);
    }
    spartan_from_pca(&index, query, query_word);
    trie_query_scratch scratch = {0};
    trie_prepare_record_lb_table(&trie, &index, query, &scratch);
    trie_query_stats stats = {0};
    const symbolic_trie_node *seed;
    file_position_type position = QUERY_RESULT_NO_POSITION;
    int ok = trie_seed_search(&index, query, query, query_word, &stats, &scratch, &seed, &position) == 0.0f &&
             seed == trie.root->children[1] && position == D && stats.checked_nodes == 1;
    free(scratch.candidates);
    N_PQUEUE = 2;
    for (maxquerythread = 1; maxquerythread <= 2; ++maxquerythread) {
        for (int finite = 0; finite < 2; ++finite) {
            const query_result result = symbolic_trie_exact_search(&index, query, query, finite ? 0.0f : FLT_MAX);
            if (result.distance != 0.0f || result.record_position != D) ok = 0;
        }
    }
    for (int d = 0; d < D; ++d) query[d] = 1.125f;
    const query_result bounded = symbolic_trie_exact_search(&index, query, query, 0.0f);
    if (bounded.distance != 0.0f || bounded.record_position != QUERY_RESULT_NO_POSITION) ok = 0;
    trie_node_destroy(trie.root);
    free(rawfile); rawfile = NULL;
    if (!ok) fprintf(stderr, "trie seed initialization lost its candidate or BSF\n");
    return ok;
}

static int check_batch_words(void) {
    enum { N = 35, FULL = 128 };
    sax_type words[N][FULL];
    float table[64][256];
    const int widths[] = {16, 17, 31, 64};
    for (int r = 0; r < N; ++r)
        for (int d = 0; d < FULL; ++d) words[r][d] = (sax_type) (r * 19 + d * 7);
    for (int d = 0; d < 64; ++d)
        for (int s = 0; s < 256; ++s) table[d][s] = (float) ((s + d * 17) % 256) * 0.000113f;
    symbolic_trie_node *leaf = trie_node_create(FULL, NULL, (trie_dimension_mask) {0, 0});
    if (!leaf) return 0;
    for (int r = 0; r < N; ++r)
        if (!trie_leaf_append(leaf, words[r], r * FULL, FULL)) return 0;
    fail_next_allocation = 1;
    int ok = trie_prepare_batch_words(leaf, 16) == 1 && leaf->batch_words == NULL;
    const int lanes = messi_record_lb_batch_lanes();
    for (size_t w = 0; w < sizeof(widths) / sizeof(*widths) && ok; ++w) {
        const int dimensions = widths[w];
        if (trie_prepare_batch_words(leaf, dimensions)) { ok = 0; break; }
        for (int r = 0; r < 48; ++r)
            for (int d = 0; d < dimensions; ++d)
                if (leaf->batch_words[((size_t) (r / 16) * dimensions + d) * 16 + r % 16] !=
                    (r < N ? words[r][d] : 0)) ok = 0;
        for (int base = 0; base < N && ok; base += lanes) {
            const float reference = messi_record_lb_table_scalar(table, words[base], dimensions, FLT_MAX);
            const float limits[] = {FLT_MAX, 0.0f, 0.25f, nextafterf(reference, -INFINITY),
                                    reference, nextafterf(reference, INFINITY)};
            for (size_t l = 0; l < sizeof(limits) / sizeof(*limits); ++l)
                for (int sparse = 0; sparse < 2; ++sparse) {
                    const int valid = N - base < lanes ? N - base : lanes;
                    const unsigned int active = ((1U << valid) - 1U) & (sparse ? 0xaaaaU : 0xffffU);
                    float distances[16];
                    const unsigned int survivors = messi_record_lb_table_batch(table,
                        leaf->batch_words + (size_t) (base / 16) * dimensions * 16 + base % 16,
                        dimensions, limits[l], active, distances);
                    unsigned int expected = 0;
                    for (int lane = 0; lane < valid; ++lane) {
                        if (!(active & (1U << lane))) continue;
                        const float sum = messi_record_lb_table_scalar(table, words[base + lane], dimensions, FLT_MAX);
                        if (sum <= limits[l]) {
                            expected |= 1U << lane;
                            if (distances[lane] != sum) ok = 0;
                        }
                    }
                    if (survivors != expected) ok = 0;
                }
        }
    }
    trie_node_destroy(leaf);
    if (!ok) fprintf(stderr, "trie transposed layout or masked bound mismatch\n");
    return ok;
}

int main(void) {
    const int dimensions[] = {1, 17, 64, 128};
    for (size_t d = 0; d < sizeof(dimensions) / sizeof(*dimensions); ++d)
        for (int split = 0; split <= 1; ++split)
            if (!check_compaction(dimensions[d], split)) {
                fprintf(stderr, "trie compaction failed: dimensions=%d split=%d\n", dimensions[d], split);
                return 1;
            }
    return check_seed_search() && check_batch_words() ? 0 : 1;
}
