#ifndef MESSI_TRIE_BATCH_H
#define MESSI_TRIE_BATCH_H

/* Small generic OpenMP scheduler used by the experimental trie leaf-batch
 * backend.  Trie-specific state remains private to trie.c. */
typedef void (*trie_batch_task_fn)(void *context, int worker, int task_index);

void trie_batch_run_tasks(int task_count, int workers,
                          trie_batch_task_fn execute, void *context);

#endif
