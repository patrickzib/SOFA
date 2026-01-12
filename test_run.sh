FILE_PATH=data_head/astro_head.bin
QUERIES_PATH=data_queries/astro_queries.bin

TS_SIZE=256
COEFF_NUMBER=64
QUERY_SIZE=10

DATASET_SIZE=10000
SAMPLE_SIZE=1000

CPU_TYPE=${82:-1}
QUEUE_NUMBER=${8:-1}

./bin/MESSI \
  --dataset $FILE_PATH \
  --dataset-size $DATASET_SIZE \
  --queries $QUERIES_PATH \
  --queries-size $QUERY_SIZE \
  --timeseries-size $TS_SIZE \
  --sample-size $SAMPLE_SIZE \
  --cpu-type $CPU_TYPE \
  --queue-number $QUEUE_NUMBER \
  --function-type 6 \
  --histogram-type 2 \
  --sfa-n-coefficients $COEFF_NUMBER \
  --is-norm
