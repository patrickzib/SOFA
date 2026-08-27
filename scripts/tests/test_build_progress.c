#include "ads/build_progress.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>

int main(void) {
    FILE *capture = tmpfile();
    if (capture == NULL) return 1;
    int saved_stderr = dup(fileno(stderr));
    if (saved_stderr < 0 || dup2(fileno(capture), fileno(stderr)) < 0) return 1;

    messi_build_progress progress;
    messi_build_progress_init(&progress);
    messi_build_progress_update(&progress, 42.9);
    messi_build_progress_update(&progress, 42.1);
    messi_build_progress_update(&progress, 99.9);
    messi_build_progress_finish(&progress);
    fflush(stderr);
    dup2(saved_stderr, fileno(stderr));
    close(saved_stderr);

    char output[1024] = {0};
    rewind(capture);
    size_t bytes = fread(output, 1, sizeof(output) - 1, capture);
    fclose(capture);
    if (bytes == 0) return 1;
    if (strstr(output, "\r>> building index: \x1b[32m0.00%\x1b[0m | elapsed") == NULL ||
        strstr(output, "\r>> building index: \x1b[32m42.90%\x1b[0m | elapsed") == NULL ||
        strstr(output, "\r>> building index: \x1b[32m99.90%\x1b[0m | elapsed") == NULL ||
        strstr(output, "\r>> building index: \x1b[32m100.00%\x1b[0m | elapsed") == NULL) return 1;

    const char *first = strstr(output, "\x1b[32m100.00%\x1b[0m | elapsed");
    return first != NULL &&
           strstr(first + 1, "\x1b[32m100.00%\x1b[0m | elapsed") == NULL ? 0 : 1;
}
