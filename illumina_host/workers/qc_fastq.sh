#!/bin/bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l jobtype=serial
#PBS -l select=1:ncpus=2:mem=5gb
#PBS -l place=pack:shared
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00

source /usr/share/Modules/init/bash
module load R

SED=/bin/sed

cd $RAW_DIR

#
# The following script runs QC on  
# a set of paired end illumina fastq
# files
#

FILE2=$(echo $FILE | $SED 's/_R1_/_R2_/')

$BIN2_DIR/SolexaQA++ analysis -d $FASTQ_DIR $FILE $FILE2
$BIN2_DIR/SolexaQA++ dynamictrim -d $FASTQ_DIR $FILE $FILE2

cd $FASTQ_DIR
$BIN2_DIR/SolexaQA++ lengthsort -d $FASTQ_DIR $FILE.trimmed $FILE2.trimmed

