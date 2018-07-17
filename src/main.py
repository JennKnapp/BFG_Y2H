import argparse
import logging
import logging.config
from create_fasta import *
from param import *
from alignment import *

if __name__ == "__main__":

    # logging #
    logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
    log = logging.getLogger("root")

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # -- create: creates fasta file from summary.csv
    # -- build: if the fasta files already exist, build index files for them
    parser.add_argument('--create', help="Summary file for making referece fasta", nargs=3)
    parser.add_argument('--build', help="Path to fasta file")
    
    # parameters for cluster
    parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("--output", help="Output path for sam files")

    args = parser.parse_args()
    
    # processing fasta file
    summary = args.create
    fasta_output = None

    if summary is not None:
        for f in summary:
            if "AD" in f:
                AD_summary = f
            elif "DB" in f:
                DB_summary = f
            else:
                fasta_output = f
        
        log.info("Creating fasta file based on summary files: %s , %s", AD_summary, DB_summary)
        log.info("Ouput fasta files will be saved into: %s", fasta_output)

        # example of create fasta for AD1DB4
        create_fasta(AD_summary, DB_summary, fasta_output, group_spec=True, AD="G1", DB="G4")
        # example of create fasta for all 
        create_fasta(AD_summary, DB_summary, fasta_output)
        log.info("Fasta files created")

    if fasta_output is None:
        fasta_output = args.build

    if fasta_output is not None:
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
        # list all samples in output dir 
        samples = os.listdir(output)
        # for each sample, do the analysis. 
        # the analysis function should take three parameters: R1, R2, group summary
        

