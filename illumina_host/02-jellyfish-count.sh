#!/bin/bash

#
# Create k-mer suffix arrays from a directory of FASTA files
#

source ./config.sh

export CWD=$PWD

PROG=`basename $0 ".sh"`
ERR_DIR="$CWD/err/$PROG"
OUT_DIR="$CWD/out/$PROG"

create_dirs "$ERR_DIR" "$OUT_DIR"

if [[ ! -d "$JELLYFISH_DIR" ]]; then
    mkdir "$JELLYFISH_DIR"
fi

cd "$REF_DIR"
COUNT=`find -maxdepth 1 -type f -name \*.fa | wc -l`
echo Found $COUNT files in \"$REF_DIR\"

i=0
for FASTA in *.fa; do
    i=$((i+1))

    export FILE=`readlink -f $FASTA`

    JOB_ID=`qsub -N jellyfish -e "$ERR_DIR/$FASTA" -o "$OUT_DIR/$FASTA" -v REF_DIR,MER_SIZE,FILE,JELLYFISH,JELLYFISH_DIR $SCRIPT_DIR/jellyfish-count.sh`

    printf "%5d: %s %s\n" $i $JOB_ID $FASTA
done

echo Submitted $i jobs for you.  Have a nice day.
