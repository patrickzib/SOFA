AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS =  ${ACLOCAL_FLAGS}
AM_CFLAGS = -fcommon
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_srcdir)/include

lib_LIBRARIES	= lib/libads.a
lib_libads_a_SOURCES = src/ads/build_progress.c src/ads/isax/isax_file_loaders.c src/ads/isax/isax_first_buffer_layer.c src/ads/isax/isax_index.c src/ads/isax/isax_node.c src/ads/isax/isax_node_buffer.c src/ads/isax/isax_node_record.c src/ads/isax/isax_node_split.c src/ads/isax/isax_query_engine.c src/ads/isax/isax_visualize_index.c src/ads/isax/inmemory_query_engine.c src/ads/isax/parallel_inmemory_query_engine.c src/ads/isax/inmemory_index_engine.c src/ads/isax/parallel_index_engine.c src/ads/isax/parallel_query_engine.c src/ads/isax/inmemory_topk_engine.c src/ads/pqueue.c src/ads/sax/sax.c src/ads/sax/ts.c src/ads/sfa/dft.c src/ads/sfa/sfa.c src/ads/sffa/frft.c src/ads/sffa/sffa.c src/ads/spartan/pca.c src/ads/spartan/spartan.c src/ads/pisa/pisa.c src/ads/calc_utils.c src/ads/trie/trie.c
lib_libads_a_CFLAGS = -I/opt/local/include $(FFTW_CFLAGS) -march=native $(SIMD_CFLAGS) -fopenmp -fcommon


bin_PROGRAMS = bin/MESSI



bin_MESSI_SOURCES = src/cli/MESSI.c
bin_MESSI_CFLAGS = -I/opt/local/include -I/usr/include/ -I/usr/local/include/ $(FFTW_CFLAGS) -march=native $(SIMD_CFLAGS) -fcommon
bin_MESSI_LDADD = lib/libads.a $(FFTW_LIBS) $(LAPACK_LIBS) $(CBLAS_LINK_LIBS) -lreadline -lm -lpthread
bin_MESSI_LDFLAGS = -L/opt/local/lib -Llib/ $(SIMD_CFLAGS) -fopenmp

check_PROGRAMS = scripts/tests/test_ed_simd scripts/tests/test_trie_partitioned_bounds scripts/tests/test_build_progress scripts/tests/test_sffa scripts/tests/test_sffa_api
TESTS = scripts/tests/test_ed_simd scripts/tests/test_trie_partitioned_bounds scripts/tests/test_build_progress scripts/tests/test_sffa scripts/tests/test_sffa_api
scripts_tests_test_ed_simd_SOURCES = scripts/tests/test_ed_simd.c
scripts_tests_test_ed_simd_LDADD = lib/libads.a -lm
scripts_tests_test_trie_partitioned_bounds_SOURCES = scripts/tests/test_trie_partitioned_bounds.c
scripts_tests_test_trie_partitioned_bounds_LDADD = lib/libads.a $(FFTW_LIBS) $(LAPACK_LIBS) $(CBLAS_LINK_LIBS) -lm
scripts_tests_test_trie_partitioned_bounds_LDFLAGS = -fopenmp
scripts_tests_test_build_progress_SOURCES = scripts/tests/test_build_progress.c
scripts_tests_test_build_progress_LDADD = lib/libads.a -lm -lpthread
scripts_tests_test_sffa_SOURCES = scripts/tests/test_sffa.c
scripts_tests_test_sffa_LDADD = lib/libads.a $(FFTW_LIBS) -lm
scripts_tests_test_sffa_api_SOURCES = scripts/tests/test_sffa_api.c src/ads/api.c
scripts_tests_test_sffa_api_LDADD = lib/libads.a $(FFTW_LIBS) $(LAPACK_LIBS) $(CBLAS_LINK_LIBS) -lm -lpthread
scripts_tests_test_sffa_api_LDFLAGS = -fopenmp

EXTRA_PROGRAMS = bench_ed_simd
bench_ed_simd_SOURCES = scripts/bench_ed_simd.c
bench_ed_simd_LDADD = lib/libads.a -lm
