#### BFG Y2H Analysis Pipeline ####

##### 1. Set up parameters in `param.py` #####

Before running the pipeline, goto `param.py` to set up variables. 

##### 2. Run this pipeline on SGE #####

`./sge.sh path_to_fastq path_to_output`

The pipeline has two parts: 

a. After you submit the job to cluster, the pipeline will first 
align the fastq files to reference and count barcodes from output sam files

b. After all the alignments finished, the pipeline will submit another batch of 
jobs for score optimizations. 

##### 3. Output #####

