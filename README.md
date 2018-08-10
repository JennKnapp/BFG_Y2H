### Before running the pipeline ###

* Please goto `param.py` to set the parameters
* Boolean parameters:
  * `ALIGNMENT` set to `True` if you want to analyze the output (see **analysis** for detail
  * `READ_COUNT` set to `True` if you want to do read counts
* `AD_summary` and `DB_summary` are used to generate gene names for read count matrix or reference files. 
  * Required col: `Group; Locus; Index; UpTag_Sequence; DnTag_Sequence`
* `REF_PATH`: reference files. All the AD reference has name `y_AD_G*`, all the DB reference has name `y_DB_G*`
### Creating fasta file and build index ###

* If you want to create fasta file from csv file, please specify a path for reference in `main.py` then run the following command:

  ` python ./src/main.py --adgroup G0 --dbgroup G0 --pfasta ./path_to_fasta/ `

  * if `adgroup` or `dbgroup` is set to `G0`, the reference will be generated for all the sequences

### bowtie alignment and read counts ###

Notes before running ..

  * **PLEASE PROVIDE ABSOLUTE PATHS**
  * The pipeline aligns AD and DB separately, please make two separate reference files. 
  * All of the fastq files must have the name: `*_R1.fastq.gz` or `*_R2.fastq.gz`

To run this on **sge** (make sure you make the output folder becore running)

  `./sge.sh fastq_file_path output_path AD_GROUP DB_GROUP`

For example:

  `./sge.sh /home/rothlab/rli/01_ngsdata/180709_merged_yAD1DB4/ /home/rothlab/rli/02_dev/08_bfg_y2h/yAD1DB4_output/ G1 G4`
   
Output

  * For each set of fastq files, one folder will be created with the sample name in the `output path`
  * The alignement log and main log will be generated in each dir
  * `*_noh.sam` is the alignment output sorted sam file without header
  * `*_rawcounts.csv` is the read count matrix
