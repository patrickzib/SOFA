#!/bin/bash
set -e

./run_scedc.sh $1 $2
./run_astro.sh $1 $2
./run_seismic.sh $1 $2
./run_sald.sh $1 $2
./run_sift1b.sh $1 $2
./run_deep1b.sh $1 $2
