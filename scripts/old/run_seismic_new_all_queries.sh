#!/bin/bash
set -e


# /vol/fob-wbib-vol2/wbi/schaefpa/messi/queries_botao

./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic.bin --in-memory --timeseries-size 256  --function-type 3 --SIMD --dataset-size 144906648 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /vol/fob-wbib-vol2/wbi/schaefpa/messi/queries_botao/noiseSeismic01.bin   --queries-size 100 --queue-number 36 --cpu-type 362 --is-norm --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000


./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic.bin --in-memory --timeseries-size 256  --function-type 3 --SIMD --dataset-size 144906648 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /vol/fob-wbib-vol2/wbi/schaefpa/messi/queries_botao/noiseSeismic001.bin   --queries-size 100 --queue-number 36 --cpu-type 362 --is-norm --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000


./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic.bin --in-memory --timeseries-size 256  --function-type 3 --SIMD --dataset-size 144906648 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /vol/fob-wbib-vol2/wbi/schaefpa/messi/queries_botao/noiseSeismic002.bin   --queries-size 100 --queue-number 36 --cpu-type 362 --is-norm --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000


./bin/MESSI --dataset /glusterfs/dfs-gfs-dist/brandjak-pub/seismic.bin --in-memory --timeseries-size 256  --function-type 3 --SIMD --dataset-size 144906648 --flush-limit 300000 --read-block 20000 --sax-cardinality 8 --queries /vol/fob-wbib-vol2/wbi/schaefpa/messi/queries_botao/noiseSeismic005.bin   --queries-size 100 --queue-number 36 --cpu-type 362 --is-norm --leaf-size 20000 --min-leaf-size 20000 --initial-lbl-size 20000



