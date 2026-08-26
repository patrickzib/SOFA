#!/usr/bin/env bash
set -Eeuo pipefail

usage() { printf 'Usage: archive_results.sh DATASET_LABEL RUN_LABEL\n'; }
die() { printf 'Error: %s\n' "$*" >&2; exit 2; }

[[ $# -eq 2 ]] || { usage >&2; exit 2; }
DATASET_LABEL=$1
RUN_LABEL=$2
LOG_ROOT=${MESSI_LOG_ROOT:-"$HOME/MESSI_logs"}
RESULTS_ROOT=${MESSI_RESULTS_ROOT:-"$HOME/MESSI_SFA_logs"}

[[ -n $DATASET_LABEL && $DATASET_LABEL != */* && $DATASET_LABEL != '.' && $DATASET_LABEL != '..' ]] || die 'invalid dataset label'
[[ -n $RUN_LABEL && $RUN_LABEL != */* && $RUN_LABEL != '.' && $RUN_LABEL != '..' ]] || die 'invalid run label'
[[ -n $RESULTS_ROOT && $RESULTS_ROOT != / ]] || die 'unsafe results root'
[[ -d $LOG_ROOT ]] || die "log root does not exist: $LOG_ROOT"

mkdir -p -- "$RESULTS_ROOT/$DATASET_LABEL"
RESULTS_ROOT_REAL=$(cd -- "$RESULTS_ROOT" && pwd -P)
DEST_PARENT_REAL=$(cd -- "$RESULTS_ROOT/$DATASET_LABEL" && pwd -P)
case "$DEST_PARENT_REAL/" in
    "$RESULTS_ROOT_REAL"/*/) ;;
    *) die 'result destination escapes the configured results root' ;;
esac

DESTINATION=$DEST_PARENT_REAL/$RUN_LABEL
rm -rf -- "$DESTINATION"
mkdir -p -- "$DESTINATION"

shopt -s nullglob dotglob
LOG_FILES=("$LOG_ROOT"/*)
if [[ ${#LOG_FILES[@]} -eq 0 ]]; then
    printf 'Warning: no logs to archive from %s\n' "$LOG_ROOT" >&2
    exit 0
fi
mv -- "${LOG_FILES[@]}" "$DESTINATION/"
