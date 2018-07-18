import argparse
import logging
import logging.config
from create_fasta import *
from param import *
from alignment import *
from suppliments import *


if __name__ == "__main__":

    # logging #
    logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
    log = logging.getLogger("root")

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # -- create: creates fasta file from summary.csv
    # -- build: if the fasta files already exist, build index files for them
    parser.add_argument("--ad", help="Summary file for AD")
    parser.add_argument("--db", help="Summary file for DB")
    parser.add_argument('--pfasta', help="Path to fasta file")
    
    # parameters for cluster
    parser.add_argument("--pfastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("--output", help="Output path for sam files")
    
    # for analysis
    parser.add_argument("--r1", help=".sam file for read one")
    parser.add_argument("--r2", help=".sam file for read two")

    args = parser.parse_args()
    
    # processing fasta file
    AD_summary = args.ad
    DB_summary = args.db
    fasta_output = args.build
    
    if MAKE_FASTA:
        
        log.info("Creating fasta file based on summary files: %s , %s", AD_summary, DB_summary)
        log.info("Ouput fasta files will be saved into: %s", fasta_output)

        # example of create fasta for AD1DB4
        create_fasta(AD_summary, DB_summary, fasta_output, group_spec=True, AD="G1", DB="G4")
        # example of create fasta for all 
        create_fasta(AD_summary, DB_summary, fasta_output)
        log.info("Fasta files created")

    if BUILD:

        list_fasta = os.listdir(fasta_output)
        log.info("Building index for: %s", ", ".join(list_fasta))

        for fasta in list_fasta:
            build_index(os.path.join(fasta_output,fasta), fasta_output)
        
        log.info("Index files built")

    ##########################################################
    ###################### Alignment #########################

    output = args.output
    
    if ALIGN:
        # input fastq is always R1
        ad = args.fastq

        # grep dir name and basename
        dir_name = os.path.dirname(ad)
        basename = os.path.basename(ad)

        # list of files in that dir
        list_fastq = os.listdir(dir_name)
        
        # get base name (sample name)
        ad_base = basename.split("_R1")[0]
        
        # find corresponding R2
        db = [i for i in list_fastq if "R2" in i and i.split("_R2")[0]==ad_base][0]
        db = os.path.join(dir_name, db)
        
        if output is None: 
            print "Please specify ouput path"
            exit(0)
        else:
            output_dir = os.path.join(output, ad_base+"/")
            os.system("mkdir -p "+output_dir)
        output = output_dir
        
        os.chdir(output)
        logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
        align_log = logging.getLogger("alignments")
        align_log.info("Starts aligning...")
        align_log.info("R1: %s", ad)
        align_log.info("R2: %s", db)
        align_log.info("Output: %s", output)
        
        bowtie_align(ad, AD_REF, output)
        bowtie_align(db, DB_REF, output)
            
        align_log.info("Alignment finished")

    # Read counts
    if ANALYSIS:
        # for each sample, do the analysis. 
        # the analysis function should take three parameters: R1, R2, group summary
        r1_sam = args.r1
        r2_sam = args.r2
        AD_summary = args.ad
        DB_summary = args.db

        logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
        analysis_log = logging.getlogger("analysis")
        
        RCmain(r1, r2, AD_genes, DB_genes)

