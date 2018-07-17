#### Before running the pipeline ####

* Please goto `param.py` to set the boolean parameters
* Boolean parameters:
  * `MAKE_FASTA` set to True if you want to make fasta reference from summary file.
  * `BUILD` set to True if you want to build index files for fasta
  * `ALIGN` set to True if you want to align the reads
  * `ANALYSIS` set to True if you want to analyze the output (see **analysis** for detail
* Other varianbles
  * 
  

#### Creating fasta file and build index ####

* If you want to create fasta file from csv file run the following command:

  ` python ./src/main.py --ad ./summary/example_AD.csv --db ./summary/example_DB.csv --pfasta ./path_to_fasta/ `

* If the fasta files already exist, you can build index by providing the fasta path:

  ` python ./src/main.py --pfasta ./path_to_fasta/ `

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

