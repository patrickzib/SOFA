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
if [ -z "$3" ]; then
    echo "Parameter dataset is empty"
    exit 1
fi
if [ -z "$4" ]; then
    echo "Parameter query-load is empty"
    exit 1
fi
if [ -z "$5" ]; then
    echo "Parameter dataset-size is empty"
    exit 1
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
FILE_PATH=/vol/tmp/schaefpa/seismic/$3
QUERIES_PATH=/vol/tmp/schaefpa/seismic/$4

TS_SIZE=256
COEFF_NUMBER=32
DATASET_SIZE=$5
SAMPLE_SIZE=1000000
QUERIES=100
LEAF_SIZE=100

COMMON_ARGS="--dataset $FILE_PATH --in-memory --timeseries-size $TS_SIZE --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block $LEAF_SIZE --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERIES --queue-number $2 --cpu-type $1 --leaf-size $LEAF_SIZE --min-leaf-size $LEAF_SIZE --initial-lbl-size $LEAF_SIZE --SIMD"
SAMPLE_ARGS="--sample-size $SAMPLE_SIZE --sample-type 3 --is-norm --tight-bound"

run_messi() {
    $MESSI_BINARY $COMMON_ARGS "$@"
}

# messi+sax+simd
run_messi --function-type 3

# messi+sfa+variance+simd+equi-depth
run_messi --function-type 4 $SAMPLE_ARGS --histogram-type 1 --coeff-number $COEFF_NUMBER

# messi+sfa+variance+simd+equi-width
run_messi --function-type 4 $SAMPLE_ARGS --histogram-type 2 --coeff-number $COEFF_NUMBER

# messi+spartan+variance+simd+equi-depth
run_messi --function-type 5 $SAMPLE_ARGS --histogram-type 1

# messi+spartan+variance+simd+equi-width
run_messi --function-type 5 $SAMPLE_ARGS --histogram-type 2
