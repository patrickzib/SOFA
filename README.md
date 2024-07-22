This is the ads isax based time-series index.

In order to install run the following commands:

Our server use 5th generation intel Xeon CPU so that the CPU flug is broadwell. If you are using other type of CPU, please change the CPU flug in Makefile.am .

./configure
make

We provide  MESSI. Note that the server that executes MESSI should support AVX2.

The MESSI-mq command (including SIMD) is 
./MESSI --dataset dataset.bin --in-memory --initial-lbl-size 2000 --leaf-size 2000 --min-leaf-size 2000 --function-type 3 --SIMD --cpu-type 82 --dataset-size 1000000 --flush-limit 300000 --read-block 20000 --queries query.bin --queries-size 10 --queue-number 2 

The MESSI-SFA command is 
./MESSI --dataset dataset.bin --in-memory --initial-lbl-size 2000 --leaf-size 2000 --min-leaf-size 2000 --function-type 4 --cpu-type 82 --dataset-size 1000000 --flush-limit 300000 --read-block 20000 --queries query.bin --queries-size 10 --queue-number 2 --sample-size 1000000 --sample-type 3 --is-norm --histogram-type 2

The MESSI-SFA (best) command is 
./MESSI --dataset dataset.bin --in-memory --initial-lbl-size 2000 --leaf-size 2000 --min-leaf-size 2000 --function-type 4 --cpu-type 82 --dataset-size 1000000 --flush-limit 300000 --read-block 20000 --queries query.bin --queries-size 10 --queue-number 2 --sample-size 1000000 --sample-type 3 --is-norm --histogram-type 2 --coeff-number 32

A help page that explains the parameters is obtained by executing ./MESSI --help

MESSI DTW is:
./MESSI --dataset dataset.bin --leaf-size 2000 --initial-lbl-size 2000 --min-leaf-size 2000 --dataset-size 10000 --flush-limit 20000000 --function-type 3 --in-memory --cpu-type 82 --queries query.bin --queries-size 10 --queue-number 24 --dtwwindowsize 12 

PS：the SIMD lower bound distance calculation code is only for 16 segments.  The real distance calculation is only for Multiples of 8.
if you use other parameters please use their SISD version (remove "_SIMD").
