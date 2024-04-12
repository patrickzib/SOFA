#!/bin/bash
set -e

# messi+sax+simd
./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_new.bin --in-memory --timeseries-size 256  --function-type 3 --SIMD --dataset-size 100000000 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_queries.bin   --queries-size 100 --queue-number $2 --cpu-type $1 --is-norm --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000


# messi+sfa+best + equi-width
./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_new.bin --in-memory --timeseries-size 256  --function-type 4 --dataset-size 100000000 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_queries.bin   --queries-size 100 --queue-number $2 --sample-size 1000000 --sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 32

# messi+sfa+best + equi-width
./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_new.bin --in-memory --timeseries-size 256  --function-type 4 --SIMD --dataset-size 100000000 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_queries.bin   --queries-size 100 --queue-number $2 --sample-size 1000000 --sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 32


# messi+sfa+best
./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_new.bin --in-memory --timeseries-size 256  --function-type 4 --dataset-size 100000000 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_queries.bin   --queries-size 100 --queue-number $2 --sample-size 1000000 --sample-type 3 --cpu-type $1 --is-norm --histogram-type 1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 32

# messi+sfa+best
./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_new.bin --in-memory --timeseries-size 256  --function-type 4 --SIMD --dataset-size 100000000 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /glusterfs/dfs-gfs-dist/brandjak-pub/seismic_queries.bin   --queries-size 100 --queue-number $2 --sample-size 1000000 --sample-type 3 --cpu-type $1 --is-norm --histogram-type 1 --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000 --coeff-number 32


