#!/bin/bash
set -e

# turingANNs - ne
./run_turingANNs_norm.sh 36 36
./copy_files.sh turingANNs_ne_4 36

# BigANN - ne
./run_bigann_norm.sh 36 36
./copy_files.sh BIGANN_ne_4 36

# spacev1B - ne
./run_spaceV1B_norm.sh 36 36
./copy_files.sh SPACEV1B_ne_4 36

# text-to-image - ne
./run_text_to_image_norm.sh 36 36
./copy_files.sh TEXTTOIMAGE_ne_4 36