#!/bin/bash
set -e

# ./run_cleanup.sh

# Define a list of items
items=(36)

# Iterate over each item in the list
for threads in "${items[@]}"
do
    echo "The current $threads is: $$threads"

    # BigANN -
    ./high_frequency/run_bigann_norm.sh $threads $threads
    
    # SALD
    ./high_frequency/run_sald.sh $threads $threads
    
    # SIFT1B
    ./high_frequency/run_sift1b.sh $threads $threads
    
    # DEEP1B
    ./high_frequency/run_deep1b.sh $threads $threads
    
    # SCEDC
    ./high_frequency/run_scedc.sh $threads $threads
    
    # ASTRO
    ./high_frequency/run_astro.sh $threads $threads
    
    # ETHC
    ./high_frequency/run_seisbench.sh $threads $threads "ETHZ.bin" "ETHZ_queries.bin" 4999932
    
    # ISC_EHB_DepthPhases
    ./high_frequency/run_seisbench.sh $threads $threads "ISC_EHB_DepthPhases.bin" "ISC_EHB_DepthPhases_queries.bin" 100000000
    
    # LenDB
    ./high_frequency/run_seisbench.sh $threads $threads "LenDB.bin" "LenDB_queries.bin" 37345260
    
    # Iquique
    ./high_frequency/run_seisbench.sh $threads $threads "Iquique.bin" "Iquique_queries.bin" 578853
    
    ## NEIC
    ./high_frequency/run_seisbench.sh $threads $threads "NEIC.bin" "NEIC_queries.bin" 93473541
    
    # OBS
    ./high_frequency/run_seisbench.sh $threads $threads "OBS.bin" "OBS_queries.bin" 15508794
    
    # OBST2024
    ./high_frequency/run_seisbench.sh $threads $threads "OBST2024.bin" "OBST2024_queries.bin" 4160286
    
    # PNW
    ./high_frequency/run_seisbench.sh $threads $threads "PNW.bin" "PNW_queries.bin" 31982766
    
    # Meier2019JGR
    ./high_frequency/run_seisbench.sh $threads $threads "Meier2019JGR.bin" "Meier2019JGR_queries.bin" 6361998
    
    # STEAD
    ./high_frequency/run_seisbench.sh $threads $threads "STEAD.bin" "STEAD_queries.bin" 87323433
    
    # TXED
    ./high_frequency/run_seisbench.sh $threads $threads "TXED.bin" "TXED_queries.bin" 35851641
    

done


