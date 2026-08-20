FILE_PATH=data_head/astro_head.bin
QUERIES_PATH=data_queries/astro_queries.bin

TS_SIZE=256
COEFF_NUMBER=64
QUERY_SIZE=100

DATASET_SIZE=10000
SAMPLE_SIZE=10000

THREADS=${1:-auto}
QUEUE_NUMBER=${2:-8}

./bin/MESSI \
  --dataset $FILE_PATH \
  --dataset-size $DATASET_SIZE \
  --queries $QUERIES_PATH \
  --queries-size $QUERY_SIZE \
  --timeseries-size $TS_SIZE \
  --sample-size $SAMPLE_SIZE \
  --threads $THREADS \
  --queue-number $QUEUE_NUMBER \
  --function-type 4 \
  --histogram-type 2 \
  --sfa-n-coefficients $COEFF_NUMBER \
  --leaf-size 1000 \
  --is-norm \
  --dynamic-root-split-variance
  # --dynamic-index 4
  # --node-split-criterion 4 \
