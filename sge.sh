#!/bin/bash

# script for submitting jobs to sge cluster
# each cluster takes two fastq files, r1 and r2

FASTQ=$1
OUTPUT=$2

# check output dir

# create the directory and run script
# attention: if outpur dir exists it will be removed
#rm -r $OUTPUT; mkdir $OUTPUT
for fastq in $FASTQ/*R1*.fastq.gz; do
  echo $fastq
  # get AD and DB group info from file name

  #/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq $fastq --output $OUTPUT
  qsub -N $(basename $fastq .fastq.gz) ./sge_sub.sh $fastq $OUTPUT
done


# check # jobs running
i=1
while [ $i -eq 1 ]; do
  COMMAND=$(qstat -u rli |wc -l)
  echo "Total jobs running... : "$COMMAND # DEBUG MEG
  if [ $COMMAND -eq 0 ]; then
    echo "Starting score optimizations ... "
    i=0
    for sample in $OUTPUT*; do
      echo $sample
      #FULLPATH=$(pwd $OUTPUT)/$sample/
      #echo $FULLPATH
      #pre=$sample"_pre"$matrix
      #med=$sample"_med"$matrix
      #high=$sample"_high"$matrix
      qsub -N $(basename $sample) ./sge_score_sub.sh $sample
      #/home/rothlab/rli/py/bin/python2.7 ./src/score.py --sample $sample 
    done
  else
    sleep 1800
  fi
done


# check # jobs running
i=1
while [ $i -eq 1 ]; do
  COMMAND=$(qstat -u rli |wc -l)
  echo "Total jobs running... : "$COMMAND # DEBUG MEG
  if [ $COMMAND -eq 0 ]; then
    i=0
    find . -mindepth 2  -name "score_optimization_all.csv" -type f -exec cat {} > $OUTPUT/summary_all.csv \;
    find . -mindepth 2  -name "score_optimization_max.csv" -type f -exec cat {} > $OUTPUT/summary_maxmcc.csv \;
  else
    sleep 600
  fi
done

echo "Job Finished!"

