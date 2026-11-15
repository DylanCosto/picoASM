#!/usr/bin/env bash
# Polishes the Autocycler paper's Flye and Autocycler assemblies with Medaka's bacterial model,
# each with the read subset it was assembled from, for the Flye+Medaka (the core of ONT's
# EPI2ME wf-bacterial-genomes) and Autocycler+Medaka rows of the baseline table.
#
# Usage: scripts/bact_medaka_paper.sh [JOBS] [THREADS_PER_JOB]
#        (scripts/bact_medaka_paper.sh --one GENOME TOOL SUBSET polishes one assembly)
#
# Reads data/bacteria/autocycler/{assemblies,reads_subsampled}; writes
# data/bacteria/bench/paper_medaka/<genome>/<tool>_medaka_<subset>.fasta, with each genome's
# reference.fasta linked beside them, so that `bact_score_paper.py` scores that folder too.
# Medaka runs in the micromamba environment $MEDAKA_ENV (default `medaka`); needs GNU parallel.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
src=$root/data/bacteria/autocycler
out=$root/data/bacteria/bench/paper_medaka

polish() {
    local genome=$1 tool=$2 subset=$3 threads=$4
    local final=$out/$genome/${tool}_medaka_$subset.fasta
    [[ -s $final ]] && return 0
    local work=$out/$genome/work_${tool}_$subset
    rm -rf "$work"
    micromamba run -n "${MEDAKA_ENV:-medaka}" medaka_consensus -i "$src/reads_subsampled/$genome/$subset.fastq" \
        -d "$src/assemblies/$genome/${tool}_$subset.fasta" -o "$work" -t "$threads" --bacteria \
        > "$out/$genome/${tool}_medaka_$subset.log" 2>&1
    mv "$work/consensus.fasta" "$final"
    rm -rf "$work"
}

if [[ ${1:-} == --one ]]; then
    polish "$2" "$3" "$4" "$5"
    exit
fi
jobs=${1:-6}
threads=${2:-4}

for genome_dir in "$src"/assemblies/*/; do
    genome=$(basename "$genome_dir")
    mkdir -p "$out/$genome"
    ln -sfn "$src/assemblies/$genome/reference.fasta" "$out/$genome/reference.fasta"
done
for genome_dir in "$src"/assemblies/*/; do
    genome=$(basename "$genome_dir")
    for tool in flye autocycler; do
        for subset in 1 2 3 4 5 6; do
            echo "$genome $tool $subset"
        done
    done
done | nice -n 19 parallel -j "$jobs" --colsep ' ' "$0" --one {1} {2} {3} "$threads"
