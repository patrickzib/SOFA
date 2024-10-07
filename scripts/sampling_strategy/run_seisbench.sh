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
if [ -z "$6" ]
   then
     echo "Samling not set using default"
     SAMPLING_FACTOR=0.01
else
     echo "Sampling factor" $6
     SAMPLING_FACTOR=$3
fi

# Configuration
MESSI_BINARY=../bin/MESSI
FILE_PATH=/vol/tmp/schaefpa/seismic/$3
QUERIES_PATH=/vol/tmp/schaefpa/seismic/$4

TS_SIZE=256
COEFF_NUMBER=32
DATASET_SIZE=$5

SAMPLE_SIZE=$(echo "$DATASET_SIZE * $SAMPLING_FACTOR" | bc)
SAMPLE_SIZE=$(echo $SAMPLE_SIZE | awk '{printf "%.0f", $0}')
QUERY_SIZE=100

# messi+sfa+variance+equi-depth
$MESSI_BINARY --dataset $FILE_PATH --in-memory --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERIES --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number $COEFF_NUMBER --SIMD

# messi+sfa+variance+simd+equi-width
$MESSI_BINARY --dataset $FILE_PATH --in-memory --timeseries-size $TS_SIZE  --function-type 4 --dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries $QUERIES_PATH --queries-size $QUERIES --queue-number $2 --sample-size $SAMPLE_SIZE --sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number $COEFF_NUMBER  --SIMD