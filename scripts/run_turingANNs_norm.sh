#!/bin/bash
set -e

if [ -z $1 ]; then
    echo "Parameter cpu-type is empty"
    exit 1
fi
if [ -z $2 ]; then
    echo "Parameter queue-number is empty"
    exit 1
fi
if [ -z "$3" ]
     then
       echo "No query-load argument supplied using default"
       QUERY="turingANNs_queries.bin"
else
       echo "Using query load" $3
       QUERY=$3
fi


# SFA parameters:
# --histogram-type
#   equi-depth splitting: 1
#   equi-width splitting: 2
# -- function-type
#   SFA: 4
#   SAX: 3
# -- sample-type
#   first-n-values sampling: 1
#   uniform sampling: 2
#   random sampling: 3
#	--sample-size
#	  Set sample size for MCB

# Configuration
MESSI_BINARY=../bin/MESSI
#FILE_PATH=/vol/tmp/schaefpa/messi_datasets/turingANNs.bin
#QUERIES_PATH=/vol/tmp/schaefpa/messi_datasets/$QUERY
FILE_PATH=/home/tmp/schaefpa/messi_datasets/turingANNs.bin
QUERIES_PATH=/home/tmp/schaefpa/messi_datasets/$QUERY

# No divisible lengths produce false results
TS_SIZE=100

COEFF_NUMBER=32
DATASET_SIZE=100000000
SAMPLE_SIZE=1000000
QUERY_SIZE=100


# messi+sax+simd
$MESSI_BINARY --dataset $FILE_PATH --apply-z-norm --is-norm --in-memory --filetype-int --timeseries-size $TS_SIZE  --function-type 3 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --cpu-type $1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --SIMD

# messi+sfa+variance+equi-depth
# $MESSI_BINARY --dataset $FILE_PATH --apply-z-norm --is-norm --in-memory --filetype-int --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 0 --SIMD

# messi+sfa+variance+simd+equi-depth
$MESSI_BINARY --dataset $FILE_PATH --apply-z-norm --is-norm --in-memory --filetype-int --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number $COEFF_NUMBER --SIMD

# messi+sfa+variance+equi-width
# $MESSI_BINARY --dataset $FILE_PATH --apply-z-norm --is-norm --in-memory --filetype-int --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 0 --SIMD

# messi+sfa+variance+simd+equi-width
$MESSI_BINARY --dataset $FILE_PATH --apply-z-norm --is-norm --in-memory --filetype-int --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number $COEFF_NUMBER  --SIMD