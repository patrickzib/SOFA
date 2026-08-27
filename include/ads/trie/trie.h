#ifndef MESSI_SYMBOLIC_TRIE_H
#define MESSI_SYMBOLIC_TRIE_H

#include "ads/isax_index.h"
#include "ads/isax_query_engine.h"

enum response symbolic_trie_build(isax_index *index, const char *path, long ts_num,
                                  int filetype_int, int apply_znorm);
query_result symbolic_trie_exact_search(isax_index *index, const ts_type *query,
                                        const ts_type *transform, float bsf);
enum response symbolic_trie_query_file(isax_index *index, const char *path, int query_count,
                                       int filetype_int, int apply_znorm, float minimum_distance);
enum response symbolic_trie_query_file_batch(isax_index *index, const char *path, int query_count,
                                             int filetype_int, int apply_znorm, float minimum_distance);
enum response symbolic_trie_query_file_leaf_batch(isax_index *index, const char *path, int query_count,
                                                  int filetype_int, int apply_znorm, float minimum_distance);
void symbolic_trie_destroy(isax_index *index);

#endif
