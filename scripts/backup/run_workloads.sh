#!/bin/bash
set -e

# SEISMIC
./run_seismic_new_workloads.sh 9 9 noiseSeismic001.bin
./copy_files.sh SEISMIC001 9
./run_seismic_new_workloads.sh 18 18 noiseSeismic001.bin
./copy_files.sh SEISMIC001 18
./run_seismic_new_workloads.sh 36 36 noiseSeismic001.bin
./copy_files.sh SEISMIC001 36

./run_seismic_new_workloads.sh 9 9 noiseSeismic002.bin
./copy_files.sh SEISMIC002 9
./run_seismic_new_workloads.sh 18 18 noiseSeismic002.bin
./copy_files.sh SEISMIC002 18
./run_seismic_new_workloads.sh 36 36 noiseSeismic002.bin
./copy_files.sh SEISMIC002 36

./run_seismic_new_workloads.sh 9 9 noiseSeismic005.bin
./copy_files.sh SEISMIC005 9
./run_seismic_new_workloads.sh 18 18 noiseSeismic005.bin
./copy_files.sh SEISMIC005 18
./run_seismic_new_workloads.sh 36 36 noiseSeismic005.bin
./copy_files.sh SEISMIC005 36

./run_seismic_new_workloads.sh 9 9 noiseSeismic01.bin
./copy_files.sh SEISMIC01 9
./run_seismic_new_workloads.sh 18 18 noiseSeismic01.bin
./copy_files.sh SEISMIC01 18
./run_seismic_new_workloads.sh 36 36 noiseSeismic01.bin
./copy_files.sh SEISMIC01 36


# SALD
./run_sald_new_workloads.sh 9 9 noiseSALD001.bin
./copy_files.sh SALD001 9
./run_sald_new_workloads.sh 18 18 noiseSALD001.bin
./copy_files.sh SALD001 18
./run_sald_new_workloads.sh 36 36 noiseSALD001.bin
./copy_files.sh SALD001 36

./run_sald_new_workloads.sh 9 9 noiseSALD002.bin
./copy_files.sh SALD002 9
./run_sald_new_workloads.sh 18 18 noiseSALD002.bin
./copy_files.sh SALD002 18
./run_sald_new_workloads.sh 36 36 noiseSALD002.bin
./copy_files.sh SALD002 36

./run_sald_new_workloads.sh 9 9 noiseSALD005.bin
./copy_files.sh SALD005 9
./run_sald_new_workloads.sh 18 18 noiseSALD005.bin
./copy_files.sh SALD005 18
./run_sald_new_workloads.sh 36 36 noiseSALD005.bin
./copy_files.sh SALD005 36

./run_sald_new_workloads.sh 9 9 noiseSALD01.bin
./copy_files.sh SALD01 9
./run_sald_new_workloads.sh 18 18 noiseSALD01.bin
./copy_files.sh SALD01 18
./run_sald_new_workloads.sh 36 36 noiseSALD01.bin
./copy_files.sh SALD01 36
