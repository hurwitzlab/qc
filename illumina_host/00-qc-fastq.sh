#!/bin/bash

source ./config.sh

export CWD="$PWD"

PROG=`basename $0 ".sh"`
ERR_DIR=$CWD/err/$PROG
OUT_DIR=$CWD/out/$PROG

create_dirs $ERR_DIR $OUT_DIR

if [[ ! -d "$FASTQ_DIR" ]]; then
    mkdir "$FASTQ_DIR"
fi

#
# QC the fastq files 
# For example:
# paired reads are in separate files:
# RNA_1_ACAGTG_L008_R1_001.fastq
# RNA_1_ACAGTG_L008_R2_001.fastq
#

cd $RAW_DIR

# send the R1 file and use the name to get the R2 file
i=0
for file in *_R1_*.fastq; do
    i=$((i+1))

    export FILE=`basename $file`

    printf '%03d: %40s' $i $FILE

    FIRST=`qsub -v SCRIPT_DIR,RAW_DIR,BIN2_DIR,BIN_DIR,FILE,FASTQ_DIR -N qc_fastq -e $ERR_DIR/$FILE -o $OUT_DIR/$FILE $SCRIPT_DIR/qc_fastq.sh`

    echo $FIRST
done
