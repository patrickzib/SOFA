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
     QUERY="bigANN_queries.bin"
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
#   SPARTAN: 5
# -- sample-type
#   first-n-values sampling: 1
#   uniform sampling: 2
#   random sampling: 3
#	--sample-size
#	  Set sample size for MCB

# Configuration
MESSI_BINARY=../bin/MESSI
FILE_PATH=/vol/tmp/schaefpa/messi_datasets/bigANN.bin
QUERIES_PATH=/vol/tmp/schaefpa/messi_datasets/$QUERY

# No divisible lengths produce false results
TS_SIZE=100
COEFF_NUMBER=32
DATASET_SIZE=100000000
SAMPLE_SIZE=1000000
QUERY_SIZE=100
LEAF_SIZE=20000

COMMON_ARGS="--dataset $FILE_PATH --apply-z-norm --filetype-int --in-memory --timeseries-size $TS_SIZE --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block $LEAF_SIZE --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --cpu-type $1 --leaf-size $LEAF_SIZE --min-leaf-size $LEAF_SIZE --initial-lbl-size $LEAF_SIZE --SIMD"
SAMPLE_ARGS="--sample-size $SAMPLE_SIZE --sample-type 3 --is-norm --tight-bound"

run_messi() {
    $MESSI_BINARY $COMMON_ARGS "$@"
}

# messi+sax+simd
run_messi --function-type 3

# messi+sfa+variance+simd+equi-depth
run_messi --function-type 4 $SAMPLE_ARGS --histogram-type 1 --sfa-n-coefficients $COEFF_NUMBER

# messi+sfa+variance+simd+equi-width
run_messi --function-type 4 $SAMPLE_ARGS --histogram-type 2 --sfa-n-coefficients $COEFF_NUMBER

# messi+spartan+variance+simd+equi-depth
run_messi --function-type 5 $SAMPLE_ARGS --histogram-type 1

# messi+spartan+variance+simd+equi-width
run_messi --function-type 5 $SAMPLE_ARGS --histogram-type 2
