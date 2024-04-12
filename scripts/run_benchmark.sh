#!/bin/bash
set -e

# BigANN
./run_bigann.sh 9 9
./copy_files.sh BIGANN 9
./run_bigann.sh 18 18
./copy_files.sh BIGANN 18
./run_bigann.sh 36 36
./copy_files.sh BIGANN 36

# spacev1B
./run_spaceV1B.sh 9 9
./copy_files.sh SPACEV1B 9
./run_spaceV1B.sh 18 18
./copy_files.sh SPACEV1B 18
./run_spaceV1B.sh 36 36
./copy_files.sh SPACEV1B 36

# text-to-image
./run_text_to_image.sh 9 9
./copy_files.sh SPACEV1B 9
./run_text_to_image.sh 18 18
./copy_files.sh SPACEV1B 18
./run_text_to_image.sh 36 36
./copy_files.sh SPACEV1B 36


# SEISMIC
#./run_seismic001.sh 9 9
#./copy_files.sh SEISMIC001 9
#./run_seismic001.sh 18 18
#./copy_files.sh SEISMIC001 18
#./run_seismic001.sh 36 36
#./copy_files.sh SEISMIC001 36

# SALD
#./run_sald001.sh 9 9
#./copy_files.sh SALD001 9
#./run_sald001.sh 18 18
#./copy_files.sh SALD001 18
#./run_sald001.sh 36 36
#./copy_files.sh SALD001 36

# SEISMIC
#./run_seismic_new.sh 9 9
#./copy_files.sh SEISMIC 9
#./run_seismic_new.sh 18 18
#./copy_files.sh SEISMIC 18
#./run_seismic_new.sh 36 36
#./copy_files.sh SEISMIC 36

# SALD
#./run_sald_new.sh 9 9
#./copy_files.sh SALD 9
#./run_sald_new.sh 18 18
#./copy_files.sh SALD 18
#./run_sald_new.sh 36 36
#./copy_files.sh SALD 36

# SIFT1B
#./run_sift1b.sh 9 9
#./copy_files.sh SIFT1b 9 
#./run_sift1b.sh 18 18
#./copy_files.sh SIFT1b 18
#./run_sift1b.sh 36 36
#./copy_files.sh SIFT1b 36

# DEEP1B
#./run_deep1b.sh 9 9
#./copy_files.sh DEEP1b 9 
#./run_deep1b.sh 18 18
#./copy_files.sh DEEP1b 18
#./run_deep1b.sh 36 36
#./copy_files.sh DEEP1b 36 

# SCEDC
#./run_scedc.sh 9 9
#./copy_files.sh SCEDC 9 
#./run_scedc.sh 18 18
#./copy_files.sh SCEDC 18
#./run_scedc.sh 36 36
#./copy_files.sh SCEDC 36 

# ASTRO
#./run_astro.sh 9 9
#./copy_files.sh ASTRO 9 
#./run_astro.sh 18 18
#./copy_files.sh ASTRO 18
#./run_astro.sh 36 36
#./copy_files.sh ASTRO 36


