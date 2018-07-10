import argparse
from create_fasta import *
from param import *
from alignment import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # if no summary file provided, skip this step
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
        
        # example of create fasta for AD1DB4
        create_fasta(AD_summary, DB_summary, fasta_output, group_spec=True, AD="G1", DB="G4")
        # example of create fasta for all 
        create_fasta(AD_summary, DB_summary, fasta_output)

    if fasta_output is None:
        fasta_output = args.build

    if fasta_output is not None:
        list_fasta = os.listdir(fasta_output)
        for fasta in list_fasta:
            build_index(os.path.join(fasta_output,fasta), fasta_output)

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
        
        main_log = os.path.join(output, ad_base+"_main.log")
        with open(main_log, "w") as log:
            log.write("R1: "+ad+"\n")
            log.write("R2: "+db+"\n")
            log.write("Output: "+output+"\n")
            log.write("Starting alignments ...\n")
        
            bowtie_align(ad, AD_REF, output)
            bowtie_align(db, DB_REF, output)
            log.write("Alignment finished ...\n")

    # Read counts
