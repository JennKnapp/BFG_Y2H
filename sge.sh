#!/bin/bash

# script for submitting jobs to sge cluster
# each cluster takes two fastq files ( r1 and r2 )

FASTQ=$1
OUTPUT=$2


for fastq in $FASTQ/*R1*.fastq.gz; do
  echo $fastq
  #/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq $fastq --output $OUTPUT 
  qsub -N $(basename $fastq .fastq.gz) ./sge_sub.sh $fastq $OUTPUT
  #break
done

