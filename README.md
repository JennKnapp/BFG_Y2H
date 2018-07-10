
#### Creating fasta file and build index ####

* If you want to create fasta file from csv file run the following command:

` python main.py --create ./summary/example_AD.csv ./summary/example_DB.csv ./path_to_fasta/ `

* If the fasta files already exist, you can build index by providing the fasta path:

` python main.py --fasta ./path_to_fasta/ `

#### bowtie alignment ####


