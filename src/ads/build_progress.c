#include "ads/build_progress.h"

#include <stdio.h>
#include <time.h>

static double messi_progress_now_seconds(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double) now.tv_sec + (double) now.tv_nsec / 1000000000.0;
}

void messi_build_progress_init(messi_build_progress *progress) {
    messi_build_progress_init_labeled(progress, "building index");
}

void messi_build_progress_init_labeled(messi_build_progress *progress, const char *label) {
    pthread_mutex_init(&progress->lock, NULL);
    progress->last_percent = -1;
    progress->label = label;
    progress->started_seconds = messi_progress_now_seconds();
    messi_build_progress_update(progress, 0.0);
}

void messi_build_progress_update(messi_build_progress *progress, double percent) {
    int whole = (int) percent;
    if (whole < 0) whole = 0;
    /* Only finish() is allowed to claim completion. */
    if (whole >= 100) whole = 99;

    pthread_mutex_lock(&progress->lock);
    if (whole > progress->last_percent) {
        const double elapsed = messi_progress_now_seconds() - progress->started_seconds;
        const double eta = percent > 0.0 ? elapsed * (100.0 - percent) / percent : -1.0;
        fprintf(stderr, "\r>> %s: \x1b[32m%.2f%%\x1b[0m | elapsed %.0fs | ETA ",
                progress->label, percent, elapsed);
        if (eta >= 0.0) fprintf(stderr, "%.0fs", eta);
        else fputs("--", stderr);
        fflush(stderr);
        progress->last_percent = whole;
    }
    pthread_mutex_unlock(&progress->lock);
}

void messi_build_progress_finish(messi_build_progress *progress) {
    pthread_mutex_lock(&progress->lock);
    const double elapsed = messi_progress_now_seconds() - progress->started_seconds;
    fprintf(stderr, "\r>> %s: \x1b[32m100.00%%\x1b[0m | elapsed %.0fs | ETA 0s\n",
            progress->label, elapsed);
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
