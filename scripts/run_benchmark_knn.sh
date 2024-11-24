#!/bin/bash
set -e

# ./run_cleanup.sh

# Define a list of items
threads=36
items=(20 50)   

# Iterate over each item in the list
for knns in "${items[@]}"
do
    echo "The current knns is: $knns"

    # BigANN -
    ./knns/run_bigann_norm_knns.sh $threads $threads $knns
    ./copy_files.sh BIGANN $knns

    # SALD
    ./knns/run_sald_knns.sh $threads $threads $knns
    ./copy_files.sh SALD $knns

    # SIFT1B
    ./knns/run_sift1b_knns.sh $threads $threads $knns
    ./copy_files.sh SIFT1b $knns

    # DEEP1B
    ./knns/run_deep1b_knns.sh $threads $threads $knns
    ./copy_files.sh DEEP1b $knns

    # SCEDC
    ./knns/run_scedc_knns.sh $threads $threads $knns
    ./copy_files.sh SCEDC $knns

    # ASTRO
    ./knns/run_astro_knns.sh $threads $threads $knns
    ./copy_files.sh ASTRO $knns

    # ETHC
    ./knns/run_seisbench_knns.sh $threads $threads "ETHZ.bin" "ETHZ_queries.bin" 4999932  $knns
    ./copy_files.sh "ETHC" $knns

    # ISC_EHB_DepthPhases
    ./knns/run_seisbench_knns.sh $threads $threads "ISC_EHB_DepthPhases.bin" "ISC_EHB_DepthPhases_queries.bin" 100000000  $knns
    ./copy_files.sh "ISC_EHB_DepthPhases" $knns

    # LenDB
    ./knns/run_seisbench_knns.sh $threads $threads "LenDB.bin" "LenDB_queries.bin" 37345260  $knns
    ./copy_files.sh "LenDB" $knns

    # Iquique
    ./knns/run_seisbench_knns.sh $threads $threads "Iquique.bin" "Iquique_queries.bin" 578853  $knns
    ./copy_files.sh "Iquique" $knns

    ## NEIC
    ./knns/run_seisbench_knns.sh $threads $threads "NEIC.bin" "NEIC_queries.bin" 93473541  $knns
    ./copy_files.sh "NEIC" $knns

    # OBS
    ./knns/run_seisbench_knns.sh $threads $threads "OBS.bin" "OBS_queries.bin" 15508794  $knns
    ./copy_files.sh "OBS" $knns

    # OBST2024
    ./knns/run_seisbench_knns.sh $threads $threads "OBST2024.bin" "OBST2024_queries.bin" 4160286  $knns
    ./copy_files.sh "OBST2024" $knns

    # PNW
    ./knns/run_seisbench_knns.sh $threads $threads "PNW.bin" "PNW_queries.bin" 31982766  $knns
    ./copy_files.sh "PNW" $knns

    # Meier2019JGR
    ./knns/run_seisbench_knns.sh $threads $threads "Meier2019JGR.bin" "Meier2019JGR_queries.bin" 6361998  $knns
    ./copy_files.sh "Meier2019JGR" $knns

    # STEAD
    ./knns/run_seisbench_knns.sh $threads $threads "STEAD.bin" "STEAD_queries.bin" 87323433  $knns
    ./copy_files.sh "STEAD" $knns

    # TXED
    ./knns/run_seisbench_knns.sh $threads $threads "TXED.bin" "TXED_queries.bin" 35851641  $knns
    ./copy_files.sh "TXED" $knns
    

done


