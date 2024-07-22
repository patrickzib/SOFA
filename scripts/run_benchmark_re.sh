s#!/bin/bash
set -e

# SEISMIC
./run_seismic.sh 36 36
./copy_files.sh SEISMIC_re 36

# SALD
./run_sald.sh 36 36
./copy_files.sh SALD_re 36

# SIFT1B
./run_sift1b.sh 36 36
./copy_files.sh SIFT1b_re 36

# DEEP1B
./run_deep1b.sh 36 36
./copy_files.sh DEEP1b_re 36

# SCEDC
./run_scedc.sh 36 36
./copy_files.sh SCEDC_re 36

# ASTRO
./run_astro.sh 36 36
./copy_files.sh ASTRO_re 36

# turingANNs - ne
./run_turingANNs_norm.sh 36 36
./copy_files.sh turingANNs_re 36

# BigANN - ne
./run_bigann_norm.sh 36 36
./copy_files.sh BIGANN_re 36

# spacev1B - ne
./run_spaceV1B_norm.sh 36 36
./copy_files.sh SPACEV1B_re 36

# text-to-image - ne
./run_text_to_image_norm.sh 36 36
./copy_files.sh TEXTTOIMAGE_re 36
