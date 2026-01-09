AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS =  ${ACLOCAL_FLAGS}
AM_CFLAGS = -fcommon
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -I$(top_srcdir)/include

lib_LIBRARIES	= lib/libads.a
lib_libads_a_SOURCES = src/ads/isax_file_loaders.c src/ads/isax_first_buffer_layer.c src/ads/isax_index.c src/ads/isax_node.c src/ads/isax_node_buffer.c src/ads/isax_node_record.c src/ads/isax_node_split.c src/ads/isax_query_engine.c src/ads/isax_visualize_index.c src/ads/pqueue.c src/ads/sax/sax.c src/ads/sax/ts.c src/ads/inmemory_query_engine.c src/ads/parallel_inmemory_query_engine.c src/ads/inmemory_index_engine.c src/ads/parallel_index_engine.c src/ads/parallel_query_engine.c src/ads/inmemory_topk_engine.c src/ads/sfa/dft.c src/ads/sfa/sfa.c src/ads/spartan/pca.c src/ads/spartan/spartan.c src/ads/calc_utils.c
lib_libads_a_CFLAGS = -I/opt/local/include $(FFTW_CFLAGS) -march=native $(SIMD_CFLAGS) -fopenmp -fcommon


bin_PROGRAMS = bin/MESSI



bin_MESSI_SOURCES = src/utils/MESSI.c
bin_MESSI_CFLAGS = -I/opt/local/include -I/usr/include/ -I/usr/local/include/ $(FFTW_CFLAGS) -march=native $(SIMD_CFLAGS) -fcommon
bin_MESSI_LDADD = lib/libads.a $(FFTW_LIBS) $(LAPACK_LIBS) -lreadline -lm -lpthread
bin_MESSI_LDFLAGS = -L/opt/local/lib -Llib/ $(SIMD_CFLAGS) -fopenmp
