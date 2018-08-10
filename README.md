### Before running the pipeline ###

* Please goto `param.py` to set the boolean parameters
* Boolean parameters:
  * `MAKE_FASTA` set to `True` if you want to make fasta reference from summary file.
  * `ALIGNMENT` set to `True` if you want to analyze the output (see **analysis** for detail
  * `READ_COUNT` set to `True` if you want to do read counts


### Creating fasta file and build index ###

* Before running the pipeline, if you want to use the function to generate reference files, please change the parameter in `param.py`
and change the 

* If you want to create fasta file from csv file, change the parameter in param.py then run the following command:

  ` python ./src/main.py --adgroup G0 --dbgroup G0 --pfasta ./path_to_fasta/ `


### bowtie alignment and read counts ###

Notes before running ..

  * **PLEASE PROVIDE ABSOLUTE PATHS**
  * The pipeline aligns AD and DB separately, please make two separate reference files
  * All of the fastq files must have the name: `*_R1.fastq.gz` or `*_R2.fastq.gz`

To run this on **sge**

  `./sgesh fastq_file_path output_path `

Output

  * For each set of fastq files, one folder will be created with the sample name in the `output path`
  * The alignement log and main log will be generated in each dir

