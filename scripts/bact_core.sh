#!/usr/bin/env bash
# Runs picoasm on the core regression set (scripts/bact_datasets_core.tsv), ONT rows then HiFi
# rows, 3 runs at once on 8 threads each, and compares each output with BASELINE's: byte-identical
# ones say so, and the scores of both labels are tabulated for the rest.
#
# Usage: scripts/bact_core.sh LABEL [BASELINE]
#        (PICOASM=path/to/picoasm picks the binary, as for bact_bench.py; freeze it first)
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
label=${1:?usage: bact_core.sh LABEL [BASELINE]}
baseline=${2:-}
list=$root/scripts/bact_datasets_core.tsv
bench=$root/data/bacteria/bench
cd "$root"
for reads in ont hifi; do
    BACT_DATASETS=$list python3 scripts/bact_bench.py --tool picoasm --label "$label" \
        --reads "$reads" --threads 8 --jobs 3
done
[ -z "$baseline" ] && exit 0
echo
changed=()
while IFS=$'\t' read -r name _; do
    case $name in '#'* | '') continue ;; esac
    if cmp -s "$bench/$name/$baseline/asm.fa" "$bench/$name/$label/asm.fa"; then
        echo "$name: identical"
    else
        echo "$name: CHANGED"
        changed+=("$name")
    fi
done < "$list"
echo "${#changed[@]} of $(grep -cv '^#' "$list") changed against $baseline"
if [ ${#changed[@]} -gt 0 ]; then
    echo
    python3 scripts/bact_compare.py <(printf '%s\n' "${changed[@]}") "$label" "$baseline"
fi
