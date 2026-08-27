#include "config.h"
#include "trie_batch.h"
#include <stddef.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void trie_batch_run_tasks(int task_count, int workers,
                          trie_batch_task_fn execute, void *context) {
    if (task_count <= 0 || execute == NULL) return;
    if (workers < 1) workers = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
    {
        const int worker = omp_get_thread_num();
#pragma omp for schedule(dynamic, 1)
        for (int task = 0; task < task_count; ++task)
            execute(context, worker, task);
    }
#else
    for (int task = 0; task < task_count; ++task)
        execute(context, 0, task);
#endif
}
