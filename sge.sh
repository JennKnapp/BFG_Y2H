#!/bin/sh

# script for submitting jobs to sge cluster
# each cluster takes two fastq files ( r1 and r2 )

FASTQ=$1
AD_REF=$2
DB_REF=$3
OUTPUT=$4

for fastq in $FATQ/*R1*.fastq.gz; do
  echo $fastq
  #qsub -V -e "" -o "" /home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq $fastq --output $OUTPUT
done

