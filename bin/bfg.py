#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import argparse
from bfg_analysis import main

parser = argparse.ArgumentParser(description='BFG-Y2H')

# make fasta from summary file (AD and DB)
# -- create: creates fasta file from summary.csv
# parser.add_argument('--pfasta', help="Path to fasta file")

# parameters for cluster
parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")
parser.add_argument("--output", help="Output path for sam files")
parser.add_argument("--mode", help="pick yeast or human or virus or hedgy", required=True)

parser.add_argument("--alignment", action="store_true", help= "turn on alignment")
parser.add_argument("--readCount", action="store_true", help= "turn on read counting")

parser.add_argument("-r", help= "redirect output files to this dir")
parser.add_argument("--cutOff", type=int, help = "assign cut off", default=20)

args = parser.parse_args()

main.main(args)
