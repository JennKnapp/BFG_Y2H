
#### Creating fasta file and build index ####

* If you want to create fasta file from csv file run the following command:

  ` python ./src/main.py --create ./summary/example_AD.csv ./summary/example_DB.csv ./path_to_fasta/ `

* If the fasta files already exist, you can build index by providing the fasta path:

  ` python ./src/main.py --fasta ./path_to_fasta/ `

#### bowtie alignment and read counts ####

* To run the pipeline with **one pair** of AD/DB fastq file, please goto `param.py` and change the 
path for AD/DB fastq file, then run the following command:

  ` python ./src/main.py --output output_path `

* To run this on **sge**

  `./sge,sh fastq_file_path output_path `
