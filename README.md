### BFG Y2H Analysis Pipeline ###

**Requirements**

* Python 3.7
* Bowtie 2 and Bowtie2 build

### Files required ###

The pipeline requires reference files and summary files before running. They can be found on GALEN: 
```
all summary files contain summary of barcode information in csv format for yeast, human and virus
path: /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_data/summary/
all reference files contain all the barcodes in fasta format
path: /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_data/reference/
```
Before running the pipeline, you need to copy everything in these two folders to your designated directory.

An example sequence in output fasta file:
```
>G1;YDL169C_BC-1;7;up
CCCTTAGAACCGAGAGTGTGGGTTAAATGGGTGAATTCAGGGATTCACTCCGTTCGTCACTCAATAA
```

### Running the pipeline  ###

0. Install from github
```
1. download the package 
2. run ./update.sh
```

1. Input arguments: 
```
usage: bfg [-h] [--fastq FASTQ] [--output OUTPUT] --mode MODE [--alignment]
           [--cutOff CUTOFF]

BFG-Y2H

optional arguments:
  -h, --help       show this help message and exit
  --fastq FASTQ    Path to all fastq files you want to analyze
  --output OUTPUT  Output path for sam files
  --mode MODE      pick yeast or human or virus or hedgy
  --alignment      turn on alignment
  --summary      path to summary files (default is set to /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_data/summary/)
  --ref      path to reference files (default is set to /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_data/reference/)
  --cutOff CUTOFF  assign cut off (default is set to 20)
```

2. All the input fastq files should have names following the format: y|hAD*DB*_GFP_(pre|med|high) (for human and yeast) 

3. Run the pipeline on GALEN
```
# this will run the pipeline using slurm         
# all the fastq files in the given folder will be processed                               
bfg --fastq /path/to/fastq_files/ --output /path/to/output_dir/  --mode yeast/human/virus/hedgy
```

### Output files  ###

a) After running the pipeline, one folder will be generated for each group pair (yAD*DB*)

b) In the output folder for each group pair, we aligned R1 and R2 separately to the reference sequences for GFP_pre, GFP_med and GFP_high.

c) `*_sorted.sam`: Raw sam files generated from bowtie2

d) `*_noh.csv`: shrinked sam files, used for scoring

e) `*_counts.csv`: barcode counts for uptags, dntags, and combined (up+dn)
