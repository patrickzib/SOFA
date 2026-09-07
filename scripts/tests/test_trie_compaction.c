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

int main(void) {
    const int dimensions[] = {1, 17, 64, 128};
    for (size_t d = 0; d < sizeof(dimensions) / sizeof(*dimensions); ++d)
        for (int split = 0; split <= 1; ++split)
            if (!check_compaction(dimensions[d], split)) {
                fprintf(stderr, "trie compaction failed: dimensions=%d split=%d\n", dimensions[d], split);
                return 1;
            }
    return 0;
}
