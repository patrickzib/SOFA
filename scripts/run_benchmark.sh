#!/bin/bash
set -e

# ./run_cleanup.sh

# Define a list of items
items=(36)  # 9 18

# Iterate over each item in the list
for threads in "${items[@]}"
do
    echo "The current $threads is: $$threads"

    # BigANN -
    ./run_bigann_norm.sh $threads $threads
    ./copy_files.sh BIGANN $threads

    # SALD
    ./run_sald.sh $threads $threads
    ./copy_files.sh SALD $threads

    # SIFT1B
    ./run_sift1b.sh $threads $threads
    ./copy_files.sh SIFT1b $threads

    # DEEP1B
    ./run_deep1b.sh $threads $threads
    ./copy_files.sh DEEP1b $threads

    # SCEDC
    ./run_scedc.sh $threads $threads
    ./copy_files.sh SCEDC $threads

    # ASTRO
    ./run_astro.sh $threads $threads
    ./copy_files.sh ASTRO $threads

    # ETHC
    ./run_seisbench.sh $threads $threads "ETHZ.bin" "ETHZ_queries.bin" 4999932
    ./copy_files.sh "ETHC" $threads

    # ISC_EHB_DepthPhases
    ./run_seisbench.sh $threads $threads "ISC_EHB_DepthPhases.bin" "ISC_EHB_DepthPhases_queries.bin" 100000000
    ./copy_files.sh "ISC_EHB_DepthPhases" $threads

    # LenDB
    ./run_seisbench.sh $threads $threads "LenDB.bin" "LenDB_queries.bin" 37345260
    ./copy_files.sh "LenDB" $threads

    # Iquique
    ./run_seisbench.sh $threads $threads "Iquique.bin" "Iquique_queries.bin" 578853
    ./copy_files.sh "Iquique" $threads

    ## NEIC
    ./run_seisbench.sh $threads $threads "NEIC.bin" "NEIC_queries.bin" 93473541
    ./copy_files.sh "NEIC" $threads

    # OBS
    ./run_seisbench.sh $threads $threads "OBS.bin" "OBS_queries.bin" 15508794
    ./copy_files.sh "OBS" $threads

    # OBST2024
    ./run_seisbench.sh $threads $threads "OBST2024.bin" "OBST2024_queries.bin" 4160286
    ./copy_files.sh "OBST2024" $threads

    # PNW
    ./run_seisbench.sh $threads $threads "PNW.bin" "PNW_queries.bin" 31982766
    ./copy_files.sh "PNW" $threads

    # Meier2019JGR
    ./run_seisbench.sh $threads $threads "Meier2019JGR.bin" "Meier2019JGR_queries.bin" 6361998
    ./copy_files.sh "Meier2019JGR" $threads

    # STEAD
    ./run_seisbench.sh $threads $threads "STEAD.bin" "STEAD_queries.bin" 87323433
    ./copy_files.sh "STEAD" $threads

    # TXED
    ./run_seisbench.sh $threads $threads "TXED.bin" "TXED_queries.bin" 35851641
    ./copy_files.sh "TXED" $threads
    

done


