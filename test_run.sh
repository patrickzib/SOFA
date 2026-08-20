FILE_PATH=data_head/astro_head.bin
QUERIES_PATH=data_queries/astro_queries.bin

TS_SIZE=256
COEFF_NUMBER=64
QUERY_SIZE=100

DATASET_SIZE=10000
SAMPLE_SIZE=10000

THREADS=${1:-auto}
QUEUE_NUMBER=${2:-8}
INDEX_TYPE=${3:-isax}

case "$INDEX_TYPE" in
  isax|trie) ;;
  *)
    printf 'Usage: %s [threads] [queue-number] [isax|trie]\n' "$0" >&2
    exit 2
    ;;
esac

args=(
  --dataset "$FILE_PATH"
  --dataset-size "$DATASET_SIZE"
  --queries "$QUERIES_PATH"
  --queries-size "$QUERY_SIZE"
  --timeseries-size "$TS_SIZE"
  --n-segments 16
  --sample-size "$SAMPLE_SIZE"
  --threads "$THREADS"
  --queue-number "$QUEUE_NUMBER"
  --function-type 4
  --histogram-type 2
  --sfa-n-coefficients "$COEFF_NUMBER"
  --leaf-size 1000
  --is-norm
  --index-type "$INDEX_TYPE"
)

if [[ $INDEX_TYPE == isax ]]; then
  args+=(--dynamic-root-split-variance)
fi

./bin/MESSI "${args[@]}"
