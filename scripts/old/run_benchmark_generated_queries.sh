#!/bin/bash
set -e

#./run_seismic.sh 36 36 "other_queries/queries_hard1p_seismic_len256_znorm.bin"
#./copy_files.sh "queries_hard1p_seismic_len256_znorm" 36

# spacev1B
#./run_spaceV1B_norm.sh 36 36
#./copy_files.sh "SPACEV1B_norm" 36

#./run_spaceV1B_norm.sh 36 36 "generated/spacev1B_noise_01.bin"
#./copy_files.sh "SPACEV1B_01" 36

./run_spaceV1B_norm.sh 36 36 "generated/spacev1B_noise_025.bin"
./copy_files.sh "SPACEV1B_ne_025" 36

./run_spaceV1B_norm.sh 36 36 "generated/spacev1B_noise_05.bin"
./copy_files.sh "SPACEV1B_ne_05" 36

./run_spaceV1B_norm.sh 36 36 "generated/spacev1B_noise_1.bin"
./copy_files.sh "SPACEV1B_ne_1" 36

#./run_spaceV1B_norm.sh 36 36 "generated/spacev1B_noise_10.bin"
#./copy_files.sh "SPACEV1B_10" 36



# text-to-image
#./run_test_to_image_norm.sh 36 36
#./copy_files.sh "TEXTTOIMAGE_norm" 36

#./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_001.bin"
#./copy_files.sh "TEXTTOIMAGE_001" 36

./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_005.bin"
./copy_files.sh "TEXTTOIMAGE_ne_005" 36

./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_01.bin"
./copy_files.sh "TEXTTOIMAGE_ne_01" 36

./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_025.bin"
./copy_files.sh "TEXTTOIMAGE_ne_025" 36

#./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_05.bin"
#./copy_files.sh "TEXTTOIMAGE_05" 36

#./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_1.bin"
#./copy_files.sh "TEXTTOIMAGE_1" 36

#./run_text_to_image_norm.sh 36 36 "generated/text-to-image_noise_10.bin"
#./copy_files.sh "TEXTTOIMAGE_10" 36


# Turing ANNs
#./run_turingANNs_norm.sh 36 36
#./copy_files.sh "turingANNs_norm" 36

#./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_001.bin"
#./copy_files.sh "turingANNs_001" 36

#./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_005.bin"
#./copy_files.sh "turingANNs_005" 36

./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_01.bin"
./copy_files.sh "turingANNs_ne_01" 36

./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_025.bin"
./copy_files.sh "turingANNs_ne_025" 36

./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_05.bin"
./copy_files.sh "turingANNs_ne_05" 36

#./run_turingANNs_norm.sh 36 36 "generated/turingANNs_noise_1.bin"
#./copy_files.sh "turingANNs_1" 36