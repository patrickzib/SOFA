#!/bin/bash
set -e

# SEISMIC
./run_seismic.sh 36 36
./copy_files.sh "SEISMIC" 36

./run_seismic.sh 36 36 "other_queries/queries_hard1p_seismic_len256_znorm.bin"
./copy_files.sh "queries_hard1p_seismic_len256_znorm" 36

./run_seismic.sh 36 36 "other_queries/queries_hard10p_seismic_len256_znorm.bin"
./copy_files.sh "queries_hard10p_seismic_len256_znorm" 36

./run_seismic.sh 36 36 "other_queries/queries_hard2p_seismic_len256_znorm.bin"
./copy_files.sh "queries_hard2p_seismic_len256_znorm" 36

./run_seismic.sh 36 36 "other_queries/queries_hard5p_seismic_len256_znorm.bin"
./copy_files.sh "queries_hard5p_seismic_len256_znorm" 36

./run_seismic.sh 36 36 "other_queries/queries_size100_seismic_len256_znorm.bin"
./copy_files.sh "queries_size100_seismic_len256_znorm" 36


# DEEP1B
#./run_deep1b.sh 36 36
#./copy_files.sh "DEEP1b" 36

#./run_deep1b.sh 36 36 "other_queries/queries_hard10p_deep1b_len96_znorm.bin"
#./copy_files.sh "queries_hard10p_deep1b_len96_znorm" 36

#./run_deep1b.sh 36 36 "other_queries/queries_hard1p_deep1b_len96_znorm.bin"
#./copy_files.sh "queries_hard1p_deep1b_len96_znorm" 36

#./run_deep1b.sh 36 36 "other_queries/queries_hard2p_deep1b_len96_znorm.bin"
#./copy_files.sh "queries_hard2p_deep1b_len96_znorm" 36

#./run_deep1b.sh 36 36 "other_queries/queries_hard5p_deep1b_len96_znorm.bin"
#./copy_files.sh "queries_hard5p_deep1b_len96_znorm" 36

#./run_deep1b.sh 36 36 "other_queries/queries_orig100_deep1b_len96_znorm.bin"
#./copy_files.sh "queries_orig100_deep1b_len96_znorm" 36



# SALD
./run_sald.sh 36 36
./copy_files.sh "SALD" 36

./run_sald.sh 36 36 "other_queries/queries_hard10p_sald_len128_znorm.bin"
./copy_files.sh "queries_hard10p_sald_len128_znorm" 36

./run_sald.sh 36 36 "other_queries/queries_hard1p_sald_len128_znorm.bin"
./copy_files.sh "queries_hard1p_sald_len128_znorm" 36

./run_sald.sh 36 36 "other_queries/queries_hard2p_sald_len128_znorm.bin"
./copy_files.sh "queries_hard2p_sald_len128_znorm" 36

./run_sald.sh 36 36 "other_queries/queries_hard5p_sald_len128_znorm.bin"
./copy_files.sh "queries_hard5p_sald_len128_znorm" 36

./run_sald.sh 36 36 "other_queries/queries_size100_sald_len128_znorm.bin"
./copy_files.sh "queries_size100_sald_len128_znorm" 36

