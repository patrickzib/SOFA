#include "ads/build_progress.h"

#include <stdio.h>

void messi_build_progress_init(messi_build_progress *progress) {
    pthread_mutex_init(&progress->lock, NULL);
    progress->last_percent = -1;
    messi_build_progress_update(progress, 0.0);
}

void messi_build_progress_update(messi_build_progress *progress, double percent) {
    int whole = (int) percent;
    if (whole < 0) whole = 0;
    /* Only finish() is allowed to claim completion. */
    if (whole >= 100) whole = 99;

    pthread_mutex_lock(&progress->lock);
    if (whole > progress->last_percent) {
        fprintf(stderr, "\r>> building index: %d%%", whole);
        fflush(stderr);
        progress->last_percent = whole;
    }
    pthread_mutex_unlock(&progress->lock);
}

void messi_build_progress_finish(messi_build_progress *progress) {
    pthread_mutex_lock(&progress->lock);
    fprintf(stderr, "\r>> building index: 100%%\n");
    fflush(stderr);
    progress->last_percent = 100;
    pthread_mutex_unlock(&progress->lock);
    pthread_mutex_destroy(&progress->lock);
}

void messi_build_progress_abort(messi_build_progress *progress) {
    pthread_mutex_lock(&progress->lock);
    if (progress->last_percent >= 0) {
        fputc('\n', stderr);
        fflush(stderr);
    }
    pthread_mutex_unlock(&progress->lock);
    pthread_mutex_destroy(&progress->lock);
}
