#!/usr/bin/env bash
DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
DATASETS=ethc,isc_ehb_depthphases,lendb,iquique,neic,obs,obst2024,pnw,meier2019jgr,stead,txed
exec "$DIR/run_suite.sh" standard --datasets "$DATASETS" "$@"
