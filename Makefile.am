AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS =  ${ACLOCAL_FLAGS}

lib_LIBRARIES	= lib/libads.a
lib_libads_a_SOURCES = src/ads/isax_file_loaders.c src/ads/isax_first_buffer_layer.c src/ads/isax_index.c src/ads/isax_node.c src/ads/isax_node_buffer.c src/ads/isax_node_record.c src/ads/isax_node_split.c src/ads/isax_query_engine.c src/ads/isax_visualize_index.c src/ads/pqueue.c src/ads/sax/sax.c src/ads/sax/ts.c src/ads/inmemory_query_engine.c src/ads/parallel_inmemory_query_engine.c src/ads/inmemory_index_engine.c src/ads/parallel_index_engine.c src/ads/parallel_query_engine.c src/ads/inmemory_topk_engine.c src/ads/dft.c src/ads/sfa.c
lib_libads_a_CFLAGS =-I/opt/local/include -Iinclude/ -I/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/include -march=native -mavx -mavx2 -msse3 -fopenmp


bin_PROGRAMS = bin/MESSI



bin_MESSI_SOURCES = src/utils/MESSI.c
bin_MESSI_CFLAGS = -I/opt/local/include -I/usr/include/ -Iinclude/ -I/usr/local/include/ -I/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/include -march=native
bin_MESSI_LDADD = -lm -lreadline -lads -lfftw3f -lm -lpthread
bin_MESSI_LDFLAGS = -L/opt/local/lib -Llib/ -L/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/.libs  -mavx -mavx2 -msse3 -fopenmp

