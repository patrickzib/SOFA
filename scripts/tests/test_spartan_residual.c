#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fftw3.h>
int query_report_interval = 0;
static int residual_fail_malloc;
static void *residual_malloc(size_t size) {
    if (residual_fail_malloc) { residual_fail_malloc = 0; return NULL; }
    return malloc(size);
}
#define malloc residual_malloc
#include "../../src/ads/trie/trie.c"
#undef malloc

static unsigned rng = 173;
static float random_value(void) {
    rng = rng*1664525U + 1013904223U;
    return (float)(rng >> 8)/16777216.0f*2.0f-1.0f;
}

static void numerical_tests(void) {
    enum { D = 64 };
    isax_index index = {0}; isax_index_settings settings = {0};
    float matrix[D*D] = {0}, mean[D], x[D], q[D], px[D], pq[D];
    double bias[D] = {0};
    settings.n_segments = settings.timeseries_size = D;
    index.settings = &settings; index.pca_dim = index.pca_components_count = D;
    index.pca_components = matrix; index.pca_mean = mean; index.pca_bias = bias;
    for (int i = 0; i < D; ++i) {
        mean[i] = random_value();
        matrix[i*D+i] = 1.0f;
        
    }
    /* Compute bias only after every mean has been initialized. */
    for (int i = 0; i < D; ++i) {
        bias[i] = 0;
        for (int j = 0; j < D; ++j) bias[i] -= (double)mean[j]*matrix[i*D+j];
    }
    spartan_residual_model model = spartan_residual_init(&index);
    assert(model.enabled);
    struct symbolic_trie_index trie = {0};
    float boundaries[D][255], *bins[D];
    settings.sax_alphabet_cardinality = 256;
    index.bins = bins; index.trie = &trie;
    trie.dimensions = D;
    for (int d = 0; d < D; ++d) {
        bins[d] = boundaries[d];
        for (int b = 0; b < 255; ++b) boundaries[d][b] = (b-127)*0.125f;
    }
    const int prefixes[] = {0, 16, 32, 63, 64};
    for (int trial = 0; trial < 2000; ++trial) {
        int k = prefixes[trial%5];
        for (int d = 0; d < D; ++d) {
            x[d] = mean[d] + random_value()*(trial%7 == 0 ? 1e-6f : 3.0f);
            q[d] = trial%11 == 0 ? x[d] : mean[d]+random_value()*3;
            if (trial%13 == 0 && d >= k) q[d] = x[d];
            if (trial%17 == 0 && d < k) q[d] = x[d];
            if (trial%19 == 0 && d >= k) q[d] = 2*mean[d]-x[d];
            if (trial%23 == 0) {
                double u = (random_value()+1.0)*0.499999+0.000001;
                x[d] = mean[d] + sqrt(-2*log(u))*cos(3.141592653589793*random_value());
            }
        }
        assert(pca_from_ts(&index, x, px) == SUCCESS);
        assert(pca_from_ts(&index, q, pq) == SUCCESS);
        spartan_residual_value vx = spartan_residual_encode(&index, &model, x, px, k);
        spartan_residual_value vq = spartan_residual_encode(&index, &model, q, pq, k);
        long double nx = 0, nq = 0, ex = 0, eq = 0, exact = 0, retained = 0;
        for (int d = 0; d < D; ++d) {
            nx += ((long double)x[d]-mean[d])*((long double)x[d]-mean[d]);
            nq += ((long double)q[d]-mean[d])*((long double)q[d]-mean[d]);
            exact += ((long double)x[d]-q[d])*((long double)x[d]-q[d]);
        }
        for (int i = 0; i < k; ++i) {
            long double zx = 0, zq = 0;
            for (int d = 0; d < D; ++d) {
                zx += matrix[i*D+d]*((long double)x[d]-mean[d]);
                zq += matrix[i*D+d]*((long double)q[d]-mean[d]);
            }
            ex += zx*zx; eq += zq*zq; retained += (zx-zq)*(zx-zq);
        }
        long double rx = sqrtl(fmaxl(0, nx-ex));
        long double rq = sqrtl(fmaxl(0, nq-eq));
        assert(fabsl(rx-vx.radius) <= 0.01L);
        assert(fabsl(rq-vq.radius) <= 0.01L);
        double gap = spartan_residual_gap(vq.radius, vx.radius, vx.radius);
        assert(retained+gap*gap <= exact+0.01L);
        trie.residual = model; trie.bound_dimensions = k;
        trie_query_scratch scratch = {0};
        trie_residual_prepare(&index, q, pq, &scratch);
        sax_type word[D]; spartan_from_pca(&index, px, word);
        float bounded = trie_residual_bound(&index, pq, word, word, word, word,
            vx.radius, vx.radius, &scratch, FLT_MAX, NULL, 0);
        assert((long double)bounded <= exact+0.01L);
    }
    q[0] = NAN;
    assert(isnan(spartan_residual_encode(&index, &model, q, pq, 16).radius));
    matrix[0] = NAN;
    assert(!spartan_residual_init(&index).enabled);
}

static void verify_layout(const symbolic_trie_node *node, const isax_index *index) {
    if (!node->leaf) {
        for (int c = 0; c < node->split_fanout; ++c) if (node->children[c]) {
            verify_layout(node->children[c], index);
        }
        return;
    }
    float projected[128];
    for (int r = 0; r < node->size; ++r) {
        const float *raw = rawfile+node->positions[r];
        assert(pca_from_ts(index, raw, projected) == SUCCESS);
        spartan_residual_value value = spartan_residual_encode(index, &index->trie->residual,
            raw, projected, index->trie->bound_dimensions);
        assert(value.radius == node->record_residuals[r]);
    }
}

static void search_tests(int capacity, int histogram, int prefix, int normalize) {
    enum { N = 4103, D = 64 };
    isax_index index = {0}; isax_index_settings settings = {0};
    float matrix[D*D] = {0}, mean[D] = {0}, boundaries[D][255], *bins[D];
    double bias[D] = {0}; sax_type cards[D];
    settings.n_segments = settings.timeseries_size = D;
    settings.function_type = 5; settings.index_type = MESSI_INDEX_TRIE;
    settings.sax_bit_cardinality = 8; settings.sax_alphabet_cardinality = 256;
    settings.max_sax_cardinalities = cards; settings.max_leaf_size = capacity;
    settings.trie_bound_dimensions = prefix; settings.trie_split_dimensions = prefix;
    settings.trie_fanout = 8; settings.trie_leaf_ivf = 4;
    settings.trie_record_mbr_suffix_bound = 1; settings.trie_residual_norm_bound = 1;
    settings.SIMD_flag = 1; settings.histogram_type = histogram;
    index.settings = &settings; index.pca_dim = index.pca_components_count = D;
    index.pca_components = matrix; index.pca_mean = mean; index.pca_bias = bias;
    index.bins = bins; index.norm_factor = 1;
    for (int d = 0; d < D; ++d) {
        matrix[d*D+d] = 1; bins[d] = boundaries[d]; cards[d] = 8;
        for (int s = 0; s < 255; ++s) boundaries[d][s] = (s-127)*0.125f;
    }
    char path[] = "/tmp/messi-residual-XXXXXX";
    int fd = mkstemp(path); assert(fd >= 0);
    FILE *file = fdopen(fd, "wb"); assert(file);
    float row[D];
    for (int r = 0; r < N; ++r) {
        for (int d = 0; d < D; ++d) row[d] = random_value()*(d < prefix ? 0.2f : 1.0f+(r%9));
        assert(fwrite(row, sizeof(float), D, file) == D);
    }
    assert(fclose(file) == 0);
    maxquerythread = 4; N_PQUEUE = 4;
    assert(symbolic_trie_build(&index, path, N, 0, normalize) == SUCCESS);
    assert(unlink(path) == 0);
    assert(index.trie->residual.enabled);
    verify_layout(index.trie->root, &index);
    for (int simd = 0; simd <= 1; ++simd) for (int streaming = 0; streaming <= 1; ++streaming)
        for (int workers = 1; workers <= 4; workers *= 4) {
            settings.SIMD_flag = simd; settings.trie_streaming_leaf_scan = streaming;
            maxquerythread = workers;
            for (int q = 0; q < 8; ++q) {
                for (int d = 0; d < D; ++d) row[d] = q == 0 ? rawfile[d] : random_value()*4;
                if (normalize && q != 0) znorm(row, D);
                float projected[D]; assert(pca_from_ts(&index, row, projected) == SUCCESS);
                double exact = INFINITY;
                for (int r = 0; r < N; ++r) {
                    double distance = 0;
                    for (int d = 0; d < D; ++d) { double delta = (double)row[d]-rawfile[r*D+d]; distance += delta*delta; }
                    exact = fmin(exact, distance);
                }
                query_result result = symbolic_trie_exact_search(&index, row, projected, FLT_MAX);
                assert(fabs(result.distance-exact) <= 1e-5*fmax(1, exact));
            }
        }
    /* Paired end-to-end searches: same queries, unchanged tree/seed ordering.
     * This is a synthetic benchmark, not a prediction for SALD. */
    float queries[16][D], projected[16][D], distances[2][16];
    for (int q = 0; q < 16; ++q) {
        for (int d = 0; d < D; ++d) queries[q][d] = random_value()*2;
        if (normalize) znorm(queries[q], D);
        assert(pca_from_ts(&index, queries[q], projected[q]) == SUCCESS);
    }
    settings.trie_streaming_leaf_scan = 1;
    for (int mode = 0; mode < 2; ++mode) {
        index.trie->residual.enabled = mode;
        unsigned long exact = 0, checks = 0, wins = 0, prunes = 0;
        double start = messi_monotonic_seconds();
        for (int q = 0; q < 16; ++q) {
            trie_query_scratch scratch = {0}; trie_query_stats stats = {0};
            sax_type word[D]; spartan_from_pca(&index, projected[q], word);
            trie_prepare_record_lb_table(index.trie, &index, projected[q], &scratch);
            const symbolic_trie_node *seed;
            float bsf = trie_seed_search(&index, queries[q], projected[q], word, &stats, &scratch, &seed, NULL);
            distances[mode][q] = trie_search_node(&index, index.trie->root, queries[q], projected[q], bsf,
                                                 &stats, seed, &scratch, NULL);
            exact += stats.exact_distances;
            checks += stats.residual_checks[2]; wins += stats.residual_wins[2]; prunes += stats.residual_prunes[2];
            if (mode) assert(distances[mode][q] == distances[0][q]);
            free(scratch.candidates);
        }
        printf("synthetic k=%d leaf=%d residual=%d queries=16 seconds=%.6f exact=%lu record_checks=%lu wins=%lu prunes=%lu\n",
            prefix, capacity, mode, messi_monotonic_seconds()-start, exact, checks, wins, prunes);
    }
    /* Failed final-layout allocation must be recoverable without touching words. */
    trie_residual_clear(index.trie->root);
    index.trie->residual_arena = calloc(N, sizeof(float)); assert(index.trie->residual_arena);
    residual_fail_malloc = 1;
    assert(!trie_residual_finish(index.trie, index.trie->root, D));
    trie_residual_clear(index.trie->root);
    free(index.trie->residual_arena); index.trie->residual_arena = NULL;
    symbolic_trie_destroy(&index); free(rawfile); rawfile = NULL;
}

int main(void) {
    numerical_tests();
    search_tests(8192, 1, 16, 0); /* IVF and partial SIMD tail. */
    search_tests(128, 2, 32, 0);  /* Multiple levels and record reordering. */
    search_tests(128, 1, 64, 0);  /* No true discarded dimensions. */
    search_tests(128, 2, 16, 1);  /* Normalization before residual encoding. */
    puts("Residual arithmetic, layout, failure cleanup and search regression tests passed");
    return 0;
}
