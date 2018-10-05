#!/bin/bash
#$ -cwd
#$ -V

/home/rothlab/rli/py/bin/python2.7 ./src/score.py --sample $1 
