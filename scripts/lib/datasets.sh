#!/usr/bin/env bash

# Dataset metadata used by run_dataset.sh. The caller supplies the profile so
# profile-specific coefficient counts remain explicit rather than duplicated in
# dozens of launch scripts.
load_dataset() {
    local dataset=${1,,}
    local profile=$2

    DATASET_ID=$dataset
    DATASET_ROOT_KIND=main
    DATASET_FILE=
    QUERY_FILE=
    TS_SIZE=
    DATASET_SIZE=
    # SFA/PISA coefficient counts include real and imaginary dimensions.
    # Prefer a 64-dimensional training pool, but respect the FFT limit of
    # timeseries_size / 2 complex coefficients (the CLI's doubled count).
    COEFF_NUMBER=64
    APPLY_Z_NORM=false
    FILETYPE_INT=false

    case "$dataset" in
        astro)
            DATASET_FILE=astro.bin; QUERY_FILE=astro_queries.bin
            TS_SIZE=256; DATASET_SIZE=100000000
            [[ $profile == high-frequency ]] && COEFF_NUMBER=128
            ;;
        bigann|bigann_norm)
            DATASET_ID=bigann
            DATASET_FILE=bigANN.bin; QUERY_FILE=bigANN_queries.bin
            TS_SIZE=100; DATASET_SIZE=100000000
            APPLY_Z_NORM=true; FILETYPE_INT=true
            [[ $profile == high-frequency ]] && COEFF_NUMBER=50
            ;;
        deep1b)
            DATASET_FILE=deep1b.bin; QUERY_FILE=deep1b_queries.bin
            TS_SIZE=96; DATASET_SIZE=100000000
            [[ $profile == high-frequency ]] && COEFF_NUMBER=48
            ;;
        sald)
            DATASET_FILE=SALD.bin; QUERY_FILE=SALD_queries.bin
            TS_SIZE=128; DATASET_SIZE=100000000
            [[ $profile == high-frequency ]] && COEFF_NUMBER=64
            ;;
        scedc)
            DATASET_FILE=SCEDC.bin; QUERY_FILE=SCEDC_queries.bin
            TS_SIZE=256; DATASET_SIZE=100000000
            if [[ $profile == standard ]]; then COEFF_NUMBER=64; fi
            if [[ $profile == high-frequency ]]; then COEFF_NUMBER=128; fi
            ;;
        sift1b)
            DATASET_FILE=sift1b.bin; QUERY_FILE=sift1b_queries.bin
            TS_SIZE=128; DATASET_SIZE=100000000
            [[ $profile == high-frequency ]] && COEFF_NUMBER=64
            ;;

        # Useful datasets migrated from scripts/old/.
        simsearchnet)
            DATASET_FILE=SimSearchNet.bin; QUERY_FILE=SimSearchNet_queries.bin
            TS_SIZE=256; DATASET_SIZE=100000000
            APPLY_Z_NORM=true; FILETYPE_INT=true
            ;;
        spacev1b)
            DATASET_FILE=spacev1B.bin; QUERY_FILE=spacev1B_queries.bin
            TS_SIZE=100; DATASET_SIZE=100000000
            APPLY_Z_NORM=true; FILETYPE_INT=true
            ;;
        text-to-image|text_to_image)
            DATASET_ID=text-to-image
            DATASET_FILE=text-to-image.bin; QUERY_FILE=text-to-image_queries.bin
            TS_SIZE=200; DATASET_SIZE=100000000
            APPLY_Z_NORM=true
            ;;
        turinganns)
            DATASET_FILE=turingANNs.bin; QUERY_FILE=turingANNs_queries.bin
            TS_SIZE=100; DATASET_SIZE=100000000
            APPLY_Z_NORM=true; FILETYPE_INT=true
            ;;
        seismic)
            DATASET_FILE=seismic.bin; QUERY_FILE=seismic_queries.bin
            TS_SIZE=256; DATASET_SIZE=100000000
            ;;

        # Named SeisBench datasets. These use a separate configurable root.
        ethc|ethz)
            DATASET_ID=ethc; DATASET_ROOT_KIND=seisbench
            DATASET_FILE=ETHZ.bin; QUERY_FILE=ETHZ_queries.bin
            TS_SIZE=256; DATASET_SIZE=4999932
            ;;
        isc_ehb_depthphases)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=ISC_EHB_DepthPhases.bin
            QUERY_FILE=ISC_EHB_DepthPhases_queries.bin
            TS_SIZE=256; DATASET_SIZE=100000000
            ;;
        lendb)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=LenDB.bin; QUERY_FILE=LenDB_queries.bin
            TS_SIZE=256; DATASET_SIZE=37345260
            ;;
        iquique)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=Iquique.bin; QUERY_FILE=Iquique_queries.bin
            TS_SIZE=256; DATASET_SIZE=578853
            ;;
        neic)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=NEIC.bin; QUERY_FILE=NEIC_queries.bin
            TS_SIZE=256; DATASET_SIZE=93473541
            ;;
        obs)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=OBS.bin; QUERY_FILE=OBS_queries.bin
            TS_SIZE=256; DATASET_SIZE=15508794
            ;;
        obst2024)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=OBST2024.bin; QUERY_FILE=OBST2024_queries.bin
            TS_SIZE=256; DATASET_SIZE=4160286
            ;;
        pnw)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=PNW.bin; QUERY_FILE=PNW_queries.bin
            TS_SIZE=256; DATASET_SIZE=31982766
            ;;
        meier2019jgr)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=Meier2019JGR.bin; QUERY_FILE=Meier2019JGR_queries.bin
            TS_SIZE=256; DATASET_SIZE=6361998
            ;;
        stead)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=STEAD.bin; QUERY_FILE=STEAD_queries.bin
            TS_SIZE=256; DATASET_SIZE=87323433
            ;;
        txed)
            DATASET_ROOT_KIND=seisbench
            DATASET_FILE=TXED.bin; QUERY_FILE=TXED_queries.bin
            TS_SIZE=256; DATASET_SIZE=35851641
            ;;
        seisbench)
            DATASET_ROOT_KIND=seisbench
            TS_SIZE=256
            ;;
        *)
            printf 'Unknown dataset: %s\n' "$1" >&2
            return 1
            ;;
    esac

    if [[ $DATASET_ROOT_KIND == seisbench && $profile == high-frequency ]]; then
        COEFF_NUMBER=128
    fi

    local max_coefficients=$((TS_SIZE / 2))
    if (( COEFF_NUMBER > max_coefficients )); then
        COEFF_NUMBER=$max_coefficients
    fi
    # The SFA/PISA coefficient count must describe complete real/imaginary
    # pairs.  All current datasets have an even maximum, but keep the
    # metadata valid for future odd-length series as well.
    if (( COEFF_NUMBER % 2 != 0 )); then
        COEFF_NUMBER=$((COEFF_NUMBER - 1))
    fi
}

active_datasets() {
    printf '%s\n' bigann sald sift1b deep1b scedc astro
}

seisbench_datasets() {
    printf '%s\n' ethc isc_ehb_depthphases lendb iquique neic obs obst2024 pnw meier2019jgr stead txed
}
