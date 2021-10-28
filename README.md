### BFG Y2H Analysis Pipeline ###

**Requirements**

* Python 3.7
* Bowtie 2 and Bowtie2 build

### Files required ###

The pipeline requires reference files before running. They can be found on GALEN: 
```
all reference files contain all the barcodes in fasta format
path: /home/rothlab/rli/02_dev/08_bfg_y2h/bfg_data/reference/
```
Before running the pipeline, you need to copy everything in these two folders to your designated directory.


#### Build new reference ###

If you need to build a new reference for your analysis, please follow:
    
1. You can refer to the create_fasta.py script to build the new fasta file 
2. Make sure the name for the sequences follows the format: `>*;ORF-BC-ID;*;up/dn`. In other words, the ORF-ID should always 
   be the second item, and the up/dn identifier should always be the last item. (see examples below)
3. Example sequences in output fasta file:
```
>G1;YDL169C_BC-1;7;up
CCCTTAGAACCGAGAGTGTGGGTTAAATGGGTGAATTCAGGGATTCACTCCGTTCGTCACTCAATAA

>G1;YMR206W_BC-1;1.0;DB;up
CCATACGAGCACATTACGGGGCTTGAGTTATATAGTCGATCCGGGCTAACTCGCATACCTCTGATAAC

>G09;56346_BC-1;24126.0;DB;dn
TCGATAGGTGCGTGTGAAGGATGTTCCCCCGGTCACCGGGCCAGTCCTCAGTCGCTCAGTCAAG
```
4. After making the fasta file, build index with bowtie2-build
`bowtie2-build filename.fasta filename`
5. Update main.py to use the summary files you generated
   * Edit parse_input_files() to add a case

### Running the pipeline  ###

* Install from pypi (recommend): `python -m pip install BFG-Y2H`

* Install and build from github, the update.sh might need to be modified before you install
```
1. download the package from github
2. inside the root folder, run ./update.sh
```

1. Input arguments: 
```
usage: bfg [-h] [--fastq FASTQ] [--output OUTPUT] --mode MODE [--alignment]
           [--ref REF] [--cutOff CUTOFF]

BFG-Y2H

optional arguments:
  -h, --help       show this help message and exit
  --fastq FASTQ    Path to all fastq files you want to analyze
  --output OUTPUT  Output path for sam files
  --mode MODE      pick yeast or human or virus or hedgy
  --alignment      turn on alignment
  --ref REF        path to all reference files
  --cutOff CUTOFF  assign cut off

```

2. All the input fastq files should have names following the format: y|hAD*DB*_GFP_(pre|med|high) (for human and yeast) 

3. Run the pipeline on GALEN
```
# this will run the pipeline using slurm         
# all the fastq files in the given folder will be processed
# run with alignment 
bfg --fastq /path/to/fastq_files/ --output /path/to/output_dir/ --mode yeast/human/virus/hedgy --alignment --ref path/to/reference

# if alignment was finished, you want to only do read counts
bfg --fastq /path/to/fastq_files/ --output /path/to/output_dir/ --mode yeast/human/virus/hedgy --ref path/to/reference
```

### Output files  ###

* After running the pipeline, one folder will be generated for each group pair (yAD*DB*)

* The folder called `GALEN_jobs` saves all the bash scripts submited to GALEN
  
* In the output folder for each group pair, we aligned R1 and R2 separately to the reference sequences for GFP_pre, GFP_med and GFP_high.

* `*_sorted.sam`: Raw sam files generated from bowtie2

* `*_noh.csv`: shrinked sam files, used for scoring

* `*_counts.csv`: barcode counts for uptags, dntags, and combined (up+dn)
