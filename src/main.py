#!/usr/bin/env python3.7

import argparse
import logging
import logging.config
import os
import glob
import re
import param
import create_fasta
import alignment
import supplements
import read_counts
import plot

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

# set global variables
###################################
# summary files are used to grep gene names, group information
# and to create fasta reference files if needed

# in the summary files, the following columns must exist: Group, Locus, Index, UpTag_Sequence, DnTag_Sequence

# summary for AD (all the genes and group)
# yeast
yAD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_AD.csv"
# human
hAD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180927_bhORFeome_AD_RL.csv"
# human with null
hvAD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180927_bhORFeome_AD_RL_withNull.csv"

# summary for DB (all the genes and group)
# yeast
yDB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_DB_AA.csv"
# human
hDB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180927_bhORFeome_DB_RL.csv"
# human with null
hvDB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180927_bhORFeome_DB_RL_withNull.csv"
# hedgy summary
heDB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20201014_hEDGY_Screen1_ORF_BC_list.csv"

# virus
vADNC = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/vADNC_withNull.csv"
vAD2u = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/vAD2u_withNull.csv"
vDBNC = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/vDBNC_withNull.csv"
vADall = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/vADall_withNull.csv"

# Directory to store all the reference sequences
hREF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/reference/h_ref/"
yREF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/reference/y_ref/"
# added for virus
hvREF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/reference/hv_ref/"
vREF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/reference/v_ref/"
# added for hedgy
heREF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/reference/h_hedgy/"

###################################

# Padding sequences used
# Padding sequences are the same for human and yeast
# DB down tags
DB_Dn1 = "TCGATAGGTGCGTGTGAAGG"
DB_Dn2 = "CCTCAGTCGCTCAGTCAAG"
# DB up tags
DB_Up1 = "CCATACGAGCACATTACGGG"
DB_Up2 = "CTAACTCGCATACCTCTGATAAC"

# AD down tags
AD_Dn1 = "CTCCAGGGTTAGGCAGATG"
AD_Dn2 = "CAATCGCACTATCCCGCTG"
# AD up tags
AD_Up1 = "CCCTTAGAACCGAGAGTGTG"
AD_Up2 = "CACTCCGTTCGTCACTCAATAA"

###################################

def check_args(arguments):

    if not os.path.isdir(arguments.o):
        os.mkdir(arguments.o)


def main(arguments):

    check_args(arguments)
    # fastq files will be provided IF and ONLY IF we are doing alignment
    # if no alignment is required and r1, r2 is provided
    # we will skip to read counts

    if arguments.r1 and arguments.r2:
        r1_csv = arguments.r1
        r2_csv = arguments.r2
        ad_base = os.path.basename(r1_csv)
        db_base = os.path.basename(r2_csv)
    else:
        print("no r1 or r2 csv files provided, proceed to alignment")
        # input fastq is always read 1
        ad = args.fastq

        # grep dir name and basename
        dir_name = os.path.dirname(ad)
        basename = os.path.basename(ad)

        # list of files in that dir
        list_fastq = os.listdir(dir_name)
        ad_base = basename.split("_R1")[0]
        # find corresponding R2
        db = [i for i in list_fastq if "R2" in i and i.split("_R2")[0]==ad_base][0]
        db = os.path.join(dir_name, db)
        db_base = db.split("_R2")[0]

        AD_GROUP, DB_GROUP, AD_REF, DB_REF = parse_input_files(arguments.mode, ad_base)

    output_dir_name = ad_base.split("_GFP_")[0]+"/"
    current_dir = os.path.dirname(os.path.abspath(__file__))
    logging.config.fileConfig(current_dir+ "logging.conf", disable_existing_loggers=False)
    log = logging.getLogger("root")

    # if output is None:
    #     if args.r is None:
    #         log.info("no output folder specified, exiting..")
    #         exit(0)
    #     else:
    #         output_dir = args.r
    #         if not os.path.isdir(output_dir):
    #             os.system("mkdir -p "+output_dir)
    # else:
    output_dir = os.path.join(arguments.o, output_dir_name)
    if not os.path.isdir(output_dir):
        os.system("mkdir -p "+output_dir)

    # output = output_dir

    if param.ALIGNMENT:

        # if sam files are provided
        if args.r1sam:
            r1_sam = args.r1sam
            r2_sam = args.r2sam
        else:

            r1_sam = alignment.bowtie_align(ad, AD_REF, output)
            r2_sam = alignment.bowtie_align(db, DB_REF, output)

        # check if sam files exist
        if not os.path.isfile(r1_sam) or not os.path.isfile(r2_sam):
            # log error
            log.error("SAM FILE DOES NOT EXIST - CHECK ALIGNMENTS")
            exit(0)

        dir_name = os.path.dirname(r1_sam)
        r1_basename = os.path.basename(r1_sam)
        r2_basename = os.path.basename(r2_sam)

        # sort r1_sam
        log.info("Sorting sam files..")
        sorted_r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_sorted.sam"))
        sort_r1 = param.SAMTOOLS+"sort -n -o "+sorted_r1+" "+r1_sam
        # sort r2_sam
        sorted_r2 = os.path.join(dir_name, r2_basename.replace(".sam", "_sorted.sam"))
        sort_r2 = param.SAMTOOLS+"sort -n -o "+sorted_r2+" "+r2_sam

        # remove headers
        r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_noh.sam"))
        r2 = os.path.join(dir_name, r2_basename.replace(".sam", "_noh.sam"))

        os.system(sort_r1)
        os.system(sort_r2)

        # remove original sam file
        os.system("rm "+r1_sam)
        os.system("rm "+r2_sam)

        log.info("Sam files are sorted")

        log.info("Removing header..")
        os.system("grep -v \"^@\" "+sorted_r1+" > "+r1)
        os.system("grep -v \"^@\" "+sorted_r2+" > "+r2)

        r1_csv = os.path.join(dir_name, r1.replace(".sam", ".csv"))
        r2_csv = os.path.join(dir_name, r2.replace(".sam", ".csv"))

        log.info("Cut file....")
        os.system("cut -f 1-5 "+r1+" > "+ r1_csv)
        os.system("cut -f 1-5 "+r2+" > "+ r2_csv)

 # remove no header sam file
        os.system("rm "+r1)
        os.system("rm "+r2)
        log.info("File shrinked")

    else:
        for f in os.listdir(output_dir):
            if ad_base in f:
                if "R1" in f and ".csv" in f:
                    r1_csv = f
                if "R2" in f and ".csv" in f:
                    r2_csv = f

    if param.READ_COUNT:

        #if not r1_csv: # if no r1_csv generated from alignment
            # user input r1 and r2 csv
        #    r1_csv = args.r1
        #    r2_csv = args.r2

        log.info("Counting reads for %s, %s", r1_csv, r2_csv)

        if args.mode == "yeast":
            AD_genes, DB_genes = supplements.read_summary(param.yAD_summary, param.yDB_summary, AD_group=AD_GROUP, DB_group=DB_GROUP)

        elif args.mode =="human":
            AD_genes, DB_genes = supplements.read_summary(param.hAD_summary, param.hDB_summary, AD_group=AD_GROUP, DB_group=DB_GROUP)

        elif args.mode == "virus":
            if "G" in AD_GROUP: # human
                AD_GROUP = AD_GROUP.split("_")[-1]
            if "G" in DB_GROUP: # human
                DB_GROUP = DB_GROUP.split("_")[-1]
            AD_genes, DB_genes = supplements.read_summary_virus(param.hvAD_summary, param.hvDB_summary, AD_group = AD_GROUP, DB_group = DB_GROUP)

        elif args.mode == "hedgy":
            AD_GROUP = AD_GROUP.split("_")[-1]
            AD_genes, DB_genes = supplements.read_summary_hedgy(param.hvAD_summary, param.heDB_summary,
                                                                AD_group=AD_GROUP)
        else:
            print("Please pride valid mode: yeast or human orvirus")
            exit(1)

        up_matrix, dn_matrix = read_counts.RCmain(r1_csv, r2_csv, AD_genes, DB_genes)
        if not args.r:
            r = "./"
        else:
            r = args.r
        if not args.cut:
            cut_off = param.cut_off
        else:
            cut_off = args.cut

        print(r)
        if not args.prefix:
            prefix = os.path.basename(r1_csv).split("R1")[0]
        else:
            prefix = args.prefix

        uptag_file = r+str(cut_off)+"_"+prefix+"_uptag_rawcounts.csv"
        dntag_file = r+str(cut_off)+"_"+prefix+"_dntag_rawcounts.csv"

        dn_matrix.to_csv(dntag_file)
        up_matrix.to_csv(uptag_file)

        combined = up_matrix + dn_matrix

        combined.to_csv(r+str(cut_off)+"_"+prefix+"_combined_counts.csv")

        basename = os.path.basename(r1_csv).split("R1")[0]
        # plot up and dn corr
        # sample_bc_corr.png
        plot.bc_corr(basename, up_matrix, dn_matrix, output =r)

        log.info("Barcode counts corr plot")



def parse_input_files(mode, ad_base):
    """
    This function is customized to read different input files and assign them different reference sequences
    Currently we can take yeast/human/virus/hedgy
    **** Please make sure the fastq filenames have the following format ****
    Input filename format:
    For yeast: yAD(1-9|M|all)DB(1-9|M|all) - 1-9 means group, M stands for Miha, all stands for all groups. e.g yAD1DBall means AD group 1 x DB all
    For human: hAD(0-9)DB(0-9), numbers stands for pooling groups. e.g hAD1DB4 means AD group 1 x DB group 4,
    For virus: h|v(AD(0-9|NC|2u|all)h|v(DB(0-9|NC|2u|all)), numbers stands for pooling groups, e.g hAD4vDBNC, human AD group 4 vs virus DBNC
    For hedgy:
    """

    if mode == "yeast":
        m = re.match(r"yAD([1-9]|M|all)DB([1-9]|M|all)", ad_base)

        AD_GROUP = "G"+m.group(1)
        DB_GROUP = "G"+m.group(2)

        AD_REF = param.yREF_PATH+"y_AD_wnull_"+AD_GROUP
        DB_REF = param.yREF_PATH+"y_DB_wnull_"+DB_GROUP

    elif mode == "human":
        m = re.match(r"hAD([0-9]+)DB([0-9]+)", ad_base)
        if int(m.group(1)) <10:
            AD_GROUP = "G0"+m.group(1)
        else:
            AD_GROUP = "G"+m.group(1)
        if int(m.group(2)) <10:
            DB_GROUP = "G0"+m.group(2)
        else:
            DB_GROUP = "G"+m.group(2)
        AD_REF = param.hREF_PATH+"h_AD_"+AD_GROUP
        DB_REF = param.hREF_PATH+"h_DB_"+DB_GROUP

    elif mode == "virus":
        # human vs virus pairwise
        m = re.match(r"(h|v)(AD([0-9]+)|ADNC|AD2u|ADall)(h|v)(DB([0-9]+)|DBNC|DB2u|DBall)", ad_base)

        if not m.group(3) is None:
            if int(m.group(3)) < 10:
                AD_GROUP = "AD_wnull_G"+m.group(3)
            else:
                AD_GROUP = "AD_wnull_G"+m.group(3)
        else:
            AD_GROUP = m.group(2)

        if not m.group(6) is None:
            if int(m.group(6)) <10:
                DB_GROUP = "DB_wnull_G"+m.group(6)
            else:
                DB_GROUP = "DB_wnull_G"+m.group(6)
        else:
            DB_GROUP = m.group(5)

        if m.group(1) == "v": #virus
            AD_REF = param.vREF_PATH+"v_"+AD_GROUP
        else:
            AD_REF = param.hvREF_PATH+"h_"+AD_GROUP

        if m.group(4) == "v": # virus
            DB_REF = param.vREF_PATH + "v_" + DB_GROUP
        else:
            DB_REF = param.hvREF_PATH + "h_" + DB_GROUP

    elif mode == "hedgy":

        m = re.match(r"hAD([0-9]+)DBhe", ad_base)
        if int(m.group(1)) <10:
            AD_GROUP = "G0"+m.group(1)
        else:
            AD_GROUP = "G"+m.group(1)

        DB_GROUP = "hedgy"

        AD_REF = param.hvREF_PATH+"h_AD_wnull_"+AD_GROUP
        DB_REF = param.heREF_PATH+"h_DB_"+DB_GROUP

    else:
        raise ValueError("Please provide valid mode: yeast, human, virus or hedgy")

    return AD_GROUP, DB_GROUP, AD_REF, DB_REF


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # -- create: creates fasta file from summary.csv
    # parser.add_argument('--pfasta', help="Path to fasta file")

    # parameters for cluster
    parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("--output", help="Output path for sam files")
    parser.add_argument("--mode", help="pick yeast or human or virus or hedgy", required=True)
    # for analysis
    parser.add_argument("--r1", help=".csv file for read one")
    parser.add_argument("--r2", help=".csv file for read two")

    parser.add_argument("--r1sam")
    parser.add_argument("--r2sam")

    parser.add_argument("--alignment", action="store_true", help= "turn on alignment")
    parser.add_argument("--read-count", action="store_true", help= "turn on read counting")

    parser.add_argument("-r", help= "redirect output files to this dir")
    parser.add_argument("-cut", type=int, help = "assign cut off")
    parser.add_argument("-prefix", type=str)
    args = parser.parse_args()
    
    ##########################################################
    ###################### Alignment #########################

    output = args.output

    # fastq files will be provided IF and ONLY IF we are doing alignment
    # if no alignment is required and r1, r2 is provided
    # we will skip to read counts
    if args.r1 and args.r2:
        r1_csv = args.r1    
        r2_csv = args.r2
        ad_base = os.path.basename(r1_csv)
        db_base = os.path.basename(r2_csv)
    else:
        print("no r1 or r2 csv files provided, proceed to alignment")
        # input fastq is always R1
        ad = args.fastq
        
        # grep dir name and basename
        dir_name = os.path.dirname(ad)
        basename = os.path.basename(ad)

        # list of files in that dir
        list_fastq = os.listdir(dir_name)
        ad_base = basename.split("_R1")[0]
        # find corresponding R2
        db = [i for i in list_fastq if "R2" in i and i.split("_R2")[0]==ad_base][0]
        db = os.path.join(dir_name, db)
        db_base = db.split("_R2")[0]
    
