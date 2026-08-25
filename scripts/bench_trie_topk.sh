#!/bin/sh
# Compare the unchanged trie 1-NN engine to the dedicated top-k engine at
# k=1.  Pass a complete bin/MESSI invocation after this script's options.
set -eu

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 bin/MESSI <normal trie query arguments>" >&2
    exit 2
fi

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/messi-trie-topk.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT HUP INT TERM

"$@" >"$work_dir/one-nn.log" 2>&1
"$@" --topk --k-size 1 >"$work_dir/topk-one.log" 2>&1

one_nn=$(sed -n 's/^>>> query wall time: \(.*\)$/\1/p' "$work_dir/one-nn.log" | tail -n 1)
topk_one=$(sed -n 's/^>>> query wall time: \(.*\)$/\1/p' "$work_dir/topk-one.log" | tail -n 1)

echo "Trie 1-NN       : ${one_nn:-not reported}"
echo "Trie top-k (k=1): ${topk_one:-not reported}"
echo "Logs were captured during the comparison; rerun the two commands directly if full output is needed."
