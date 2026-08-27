#ifndef MESSI_BUILD_PROGRESS_H
#define MESSI_BUILD_PROGRESS_H

#include <pthread.h>

/* A small, process-local terminal reporter for long-running index builders.
 * Callers may update it concurrently; it writes only when a whole percentage
 * is reached and reserves 100% for messi_build_progress_finish(). */
typedef struct messi_build_progress {
    pthread_mutex_t lock;
    int last_percent;
} messi_build_progress;

void messi_build_progress_init(messi_build_progress *progress);
void messi_build_progress_update(messi_build_progress *progress, double percent);
void messi_build_progress_finish(messi_build_progress *progress);
void messi_build_progress_abort(messi_build_progress *progress);

#endif
