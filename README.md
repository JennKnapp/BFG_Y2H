### BFG Y2H Analysis Pipeline ###

**Requirements**

* Python 2.7
* Bowtie 2 and Bowtie2 build (change path in `param.py`)
* Samtools1.4 (change path in `param.py`)

#### 1. Set up parameters in `param.py` ####

Before running the pipeline, goto `param.py` to set up variables. 

#### 2. Run this pipeline on SGE ####

More about [Sun Grid Engine](http://gridscheduler.sourceforge.net/howto/GridEngineHowto.html) 

You probably want to change the path to python2.7 in `sge_sub.sh` and `sge_score_sub.sh`

`./sge.sh -f path_to_fastq -o path_to_output`

**Note:** if only `-o` is provided, the pipeline will only do the score optimizations. 

To run score optimization on single group:

`python ./src/score.py --sample path_to_sample_output_folder`

The pipeline has two parts: 

a. After you submit the job to cluster, the pipeline will first 
align the fastq files to reference and count barcodes from output sam files

b. After all the alignments finished, the pipeline will submit another batch of 
jobs for score optimizations. 

#### 3. Output ####

For `2.a.` A folder will be generated for each group. Including all the alignment summary files and raw barcode
counts

For `2.b.` Summary of top 5 mcc for each group can be found in the output folder (`summary_maxmcc.csv`)

#### 4. Special case ####

If you want to run the code on just one sample:

`/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq fastq_read_one --output output_dir`
