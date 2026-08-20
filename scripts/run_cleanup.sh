#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd -P)
LOG_ROOT=${MESSI_LOG_ROOT:-"$HOME/MESSI_logs"}
INDEX_ROOT=${MESSI_INDEX_ROOT:-"$HOME/index"}

clear_directory() {
    local directory=$1 label=$2
    [[ -n $directory && $directory != / ]] || { printf 'Refusing unsafe %s path: %s\n' "$label" "$directory" >&2; exit 2; }
    mkdir -p -- "$directory"
    local canonical
    canonical=$(cd -- "$directory" && pwd -P)
    [[ $canonical != / ]] || { printf 'Refusing to clear /\n' >&2; exit 2; }
    [[ $canonical != "$HOME" ]] || { printf 'Refusing to clear the home directory\n' >&2; exit 2; }
    [[ $canonical != "$REPO_ROOT" && $canonical != "$SCRIPT_DIR" ]] || {
        printf 'Refusing to clear the repository or scripts directory\n' >&2
        exit 2
    }
    find "$canonical" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
}

clear_directory "$LOG_ROOT" 'log root'
clear_directory "$INDEX_ROOT" 'index root'
