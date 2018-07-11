
#### Creating fasta file and build index ####

* If you want to create fasta file from csv file run the following command:

  ` python ./src/main.py --create ./summary/example_AD.csv ./summary/example_DB.csv ./path_to_fasta/ `

* If the fasta files already exist, you can build index by providing the fasta path:

  ` python ./src/main.py --fasta ./path_to_fasta/ `

#### bowtie alignment and read counts ####


Notes before running ..

  * **PLEASE PROVIDE ABSOLUTE PATHS**
  * The pipeline aligns AD and DB separately, please make two separate reference files
  * All of the fastq files must have the name: `*_R1.fastq.gz` or `*_R2.fastq.gz`


To run the pipeline with **one pair** of AD/DB fastq file, please goto `param.py` and change the 
path for AD/DB fastq file, then run the following command:

  ` python ./src/main.py --output output_path `

To run this on **sge**

  `./sgesh fastq_file_path output_path `

Output

  * For each set of fastq files, one folder will be created with the sample name in the `output path`
  * The alignement log and main log will be generated in each dir

