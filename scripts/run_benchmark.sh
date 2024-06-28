#!/bin/bash
set -e

# ./run_cleanup.sh

# Define a list of items
items=(9 18 36)

# Iterate over each item in the list
for threads in "${items[@]}"
do
    echo "The current $threads is: $$threads"

    # SimSearchNet -
    #./run_SimSearchNet_norm.sh $threads $threads
    #./copy_files.sh SimSearchNet $threads

    # turingANNs -
    #./run_turingANNs_norm.sh 3$threads $threads
    #./copy_files.sh turingANNs $threads

    # text-to-image -
    #./run_text_to_image_norm.sh $threads $threads
    #./copy_files.sh TEXTTOIMAGE $threads

    # SEISMIC
    #./run_seismic.sh $threads $threads
    #./copy_files.sh SEISMIC $threads

    # spacev1B -
    ##./run_spaceV1B_norm.sh $threads $threads
    ##./copy_files.sh SPACEV1B $threads


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

    # GEOFON
    #./run_seisbench.sh $threads $threads "GEOFON.bin" "GEOFON_queries.bin" 38080393
    #./copy_files.sh "GEOFON" $threads

    # InstanceCounts
    #./run_seisbench.sh $threads $threads "InstanceCounts.bin" "InstanceCounts_queries.bin" 100000000
    #./copy_files.sh "InstanceCounts" $threads


    # LFEStacksSanAndreasShelly2017
    #./run_seisbench.sh $threads $threads "LFEStacksSanAndreasShelly2017.bin" "LFEStacksSanAndreasShelly2017_queries.bin" 215511
    #./copy_files.sh "LFEStacksSanAndreasShelly2017" $threads

    # LFEStacksCascadiaBostock2015
    #./run_seisbench.sh $threads $threads "LFEStacksCascadiaBostock2015.bin" "LFEStacksCascadiaBostock2015_queries.bin" 169890
    #./copy_files.sh "LFEStacksCascadiaBostock2015" $threads

    # LFEStacksMexicoFrank2014
    #./run_seisbench.sh $threads $threads "LFEStacksMexicoFrank2014.bin" "LFEStacksMexicoFrank2014_queries.bin" 1111683
    #./copy_files.sh "LFEStacksMexicoFrank2014" $threads

    # Ross2018GPD
    #./run_seisbench.sh $threads $threads "Ross2018GPD.bin" "Ross2018GPD_queries.bin" 14320950
    #./copy_files.sh "Ross2018GPD" $threads

    # Ross2018JGRFM
    #./run_seisbench.sh $threads $threads "Ross2018JGRFM.bin" "Ross2018JGRFM_queries.bin" 29082888
    #./copy_files.sh "Ross2018JGRFM" $threads

    # Ross2018JGRPick
    #./run_seisbench.sh $threads $threads "Ross2018JGRPick.bin" "Ross2018JGRPick_queries.bin" 29082888
    #./copy_files.sh "Ross2018JGRPick" $threads

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


