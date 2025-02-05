#!/bin/bash
set -e

# ./run_cleanup.sh

# Define a list of items
threads=36
items=(0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5)

# Iterate over each item in the list
for sampling_factor in "${items[@]}"
do
    echo "The current threads is: $threads"
    echo "The current factor is: $sampling_factor"

    # BigANN -
    ./sampling_strategy/run_bigann_norm.sh $threads $threads $sampling_factor
    ./copy_files.sh BIGANN $sampling_factor

    # SALD
    ./sampling_strategy/run_sald.sh $threads $threads $sampling_factor
    ./copy_files.sh SALD $sampling_factor

    # SIFT1B
    ./sampling_strategy/run_sift1b.sh $threads $threads $sampling_factor
    ./copy_files.sh SIFT1b $sampling_factor

    # DEEP1B
    ./sampling_strategy/run_deep1b.sh $threads $threads $sampling_factor
    ./copy_files.sh DEEP1b $sampling_factor

    # SCEDC
    ./sampling_strategy/run_scedc.sh $threads $threads $sampling_factor
    ./copy_files.sh SCEDC $sampling_factor

    # ASTRO
    ./sampling_strategy/run_astro.sh $threads $threads $sampling_factor
    ./copy_files.sh ASTRO $sampling_factor

    # ETHC
    ./sampling_strategy/run_seisbench.sh $threads $threads "ETHZ.bin" "ETHZ_queries.bin" 4999932
    ./copy_files.sh "ETHC" $sampling_factor

    # ISC_EHB_DepthPhases
    ./sampling_strategy/run_seisbench.sh $threads $threads "ISC_EHB_DepthPhases.bin" "ISC_EHB_DepthPhases_queries.bin" 100000000
    ./copy_files.sh "ISC_EHB_DepthPhases" $sampling_factor

    # LenDB
    ./sampling_strategy/run_seisbench.sh $threads $threads "LenDB.bin" "LenDB_queries.bin" 37345260
    ./copy_files.sh "LenDB" $sampling_factor

    # Iquique
    ./sampling_strategy/run_seisbench.sh $threads $threads "Iquique.bin" "Iquique_queries.bin" 578853
    ./copy_files.sh "Iquique" $sampling_factor

    ## NEIC
    ./sampling_strategy/run_seisbench.sh $threads $threads "NEIC.bin" "NEIC_queries.bin" 93473541
    ./copy_files.sh "NEIC" $sampling_factor

    # OBS
    ./sampling_strategy/run_seisbench.sh $threads $threads "OBS.bin" "OBS_queries.bin" 15508794
    ./copy_files.sh "OBS" $sampling_factor

    # OBST2024
    ./sampling_strategy/run_seisbench.sh $threads $threads "OBST2024.bin" "OBST2024_queries.bin" 4160286
    ./copy_files.sh "OBST2024" $sampling_factor

    # PNW
    ./sampling_strategy/run_seisbench.sh $threads $threads "PNW.bin" "PNW_queries.bin" 31982766
    ./copy_files.sh "PNW" $sampling_factor

    # Meier2019JGR
    ./sampling_strategy/run_seisbench.sh $threads $threads "Meier2019JGR.bin" "Meier2019JGR_queries.bin" 6361998
    ./copy_files.sh "Meier2019JGR" $sampling_factor

    # STEAD
    ./sampling_strategy/run_seisbench.sh $threads $threads "STEAD.bin" "STEAD_queries.bin" 87323433
    ./copy_files.sh "STEAD" $sampling_factor

    # TXED
    ./sampling_strategy/run_seisbench.sh $threads $threads "TXED.bin" "TXED_queries.bin" 35851641
    ./copy_files.sh "TXED" $sampling_factor
    

done


