
#### Creating fasta file and build index ####

* If you want to create fasta file from csv file run the following command:

  ` python ./src/main.py --create ./summary/example_AD.csv ./summary/example_DB.csv ./path_to_fasta/ `

* If the fasta files already exist, you can build index by providing the fasta path:

  ` python ./src/main.py --fasta ./path_to_fasta/ `

#### bowtie alignment and read counts ####

* Set the parameters in `param.py`

* To run the pipeline with **one pair** of AD/DB fastq file, please run the following command:

  ` python ./src/main.py --ad ad.fastq ad_ref --db db.fastq db_ref --output output_path `

* To run this on **sge**

  `./sge,sh fastq_file_path ad_ref db_ref output_path `
