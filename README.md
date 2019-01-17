### BFG Y2H Analysis Pipeline ###

**Requirements**

* Python 2.7
* Bowtie 2 and Bowtie2 build (change path in `param.py`)
* Samtools1.4 (change path in `param.py`)

#### Create fasta files from summary table ####

a. In param.py, set full path to `AD/DB_summary.csv`

b. In param.py, set path to reference (`REF_PATH`:A folder to save the output fasta filees)

c. In param.py, set padding sequences (optional)

d. In the summary table (AD and DB), the following columns are required: Group, Locus, Index

e. Run the command: `python create_fasta.py`

f. An example sequence in output fasta file:
```
>G1;YDL169C_BC-1;7;up
CCCTTAGAACCGAGAGTGTGGGTTAAATGGGTGAATTCAGGGATTCACTCCGTTCGTCACTCAATAA
```

#### Running the pipeline  ####

a. Set parameters in param.py (refer to the comments)

b. Input argument: 
```
--fastq /path/to/fastq_file    # we take R1 as input file, there should be an _R2 file in the folder 
--output /path/to/output_dir/  # make this directory before you run the pipeline
```

c. To run the pipeline on one pair of fastq files
`python ./src/main.py --fastq /path/to/fastq_file --output /path/to/output_dir/`

(NOTE: we assume all the reference files are named in the format: y_AD* or y_DB*)

(NOTE: we assume all the fastq files have "yAD*DB*" in the file name)

(NOTE: we use these filenames to match group with reference group)

d. To run the pipeline on sge
```
# this will run the pipeline using sun grid engine                                        
# all the fastq files in the given folder will be processed                               
# The script checks how many jobs are running for this user, change user name accordingly 
./sge.sh -f /path/to/fastq_files/ -o /path/to/output_dir/                                 
```
(NOTE: More about [Sun Grid Engine](http://gridscheduler.sourceforge.net/howto/GridEngineHowto.html)) 

#### 4. Special case ####

If you want to run the code on just one sample:

`/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq fastq_read_one --output output_dir`
