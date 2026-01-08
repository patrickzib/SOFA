FILE_PATH=data_head/astro_head.bin
QUERIES_PATH=data_queries/astro_queries.bin
TS_SIZE=256

COEFF_NUMBER=16
DATASET_SIZE=1000
SAMPLE_SIZE=1000
QUERY_SIZE=10
LEAF_SIZE=200
MIN_LEAF_SIZE=200
INITIAL_LBL_SIZE=$LEAF_SIZE

CPU_TYPE=${82:-1}
QUEUE_NUMBER=${8:-1}

./bin/MESSI \
  --dataset $FILE_PATH \
  --in-memory \
  --timeseries-size $TS_SIZE \
  --function-type 5 \
  --dataset-size $DATASET_SIZE \
  --flush-limit 300000 \
  --read-block 200 \
  --sax-cardinality 8 \
  --queries $QUERIES_PATH \
  --queries-size $QUERY_SIZE \
  --queue-number $QUEUE_NUMBER \
  --sample-size $SAMPLE_SIZE \
  --sample-type 3 \
  --cpu-type $CPU_TYPE \
  --is-norm \
  --histogram-type 2 \
  --leaf-size $LEAF_SIZE \
  --min-leaf-size $MIN_LEAF_SIZE \
  --initial-lbl-size $INITIAL_LBL_SIZE \
  --coeff-number $COEFF_NUMBER
