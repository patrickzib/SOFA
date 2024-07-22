# Connect to HU-Servers via remote
xfreerdp /u:schaefpa /v:gruenau1.informatik.hu-berlin.de

# Datasets
https://ln5.sync.com/dl/0b8135230/39vxx8su-tkfi7t2s-dgsvh8rp-k8ixcs8p/view/text/17004548360004?sync_id=16907535250004


# FFTW:
# /vol/home-vol3/wbi/schaefpa/fftw-3.3.10

# ./configure --enable-threads --enable-shared --enable-float

./configure --enable-float --enable-threads --enable-sse && make -j16
#make
#make install


# MESSI:

# edit Makefile.am
# Modify to add fftw3
bin_MESSI_CFLAGS = -I/opt/local/include -Iinclude/ -I/usr/local/include/ -I/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/include -march=native 

lib_libads_a_CFLAGS =-I/opt/local/include -Iinclude/ -I/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/include -march=native -mavx -mavx2 -msse3 -fopenmp

bin_MESSI_LDFLAGS = -L/opt/local/lib -Llib/ -L/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/.lib -mavx -mavx2 -msse3 -fopenmp


# Building
export LDFLAGS=-L/vol/home-vol3/wbi/schaefpa/fftw-3.3.10/.lib
./configure
make clean
make



