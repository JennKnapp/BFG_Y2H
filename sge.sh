#!/bin/bash

# script for submitting jobs to sge cluster
# each cluster takes two fastq files, r1 and r2

FASTQ=$1
OUTPUT=$2
AD=$3
DB=$4

# check output dir


# create the directory and run script
# attention: if outpur dir exists it will be removed
#rm -r $OUTPUT; mkdir $OUTPUT
for fastq in $FASTQ/*R1*.fastq.gz; do
  echo $fastq
  #/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq $fastq --output $OUTPUT
  qsub -N $(basename $fastq .fastq.gz) ./sge_sub.sh $fastq $OUTPUT $AD $DB
done


