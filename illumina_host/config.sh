#!/bin/bash

export RAW_DIR="/rsgrps/bhurwitz/hurwitzlab/data/raw/Doetschman_20111007/L008/Project_RNA_1/Sample_RNA_1"
export BASE_DIR="/rsgrps/bhurwitz/hurwitzlab/data/processed/Doetschman_20111007/L008/Project_RNA_1/Sample_RNA_1"
export FASTQ_DIR="$BASE_DIR/data/fastq"
export FASTA_DIR="$BASE_DIR/data/fasta"
export REF_DIR="/rsgrps/bhurwitz/hurwitzlab/data/reference/mouse_genome/20141111"
export SUFFIX_DIR="$BASE_DIR/data/suffix"
export KMER_DIR="$BASE_DIR/data/kmer"
export JELLYFISH_DIR="$BASE_DIR/data/jellyfish"
export SCRIPT_DIR="$BASE_DIR/scripts/workers"
export COUNT_DIR="$BASE_DIR/data/counts"
export MER_SIZE=20
export BIN_DIR="/rsgrps/bhurwitz/bin"
export JELLYFISH="$BIN_DIR/jellyfish"

function create_dirs {
    for dir in $1 $2; do
        if [ -d "$dir" ]; then
            rm -rf "$dir/*"
        else
            mkdir -p "$dir"
        fi
    done
}
