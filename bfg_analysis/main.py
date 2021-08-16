#!/usr/bin/env python3.7

import argparse
import logging
import logging.config
import os
import glob
import re
from legacy import param
import alignment
import read_counts

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

    if not os.path.isdir(arguments.fastq):
        raise NotADirectoryError("Cannot find dir: {arguments.fastq}")


def main(arguments):

    check_args(arguments)
    # fastq files will be provided IF and ONLY IF we are doing alignment
    # if no alignment is required and r1, r2 is provided
    # we will skip to read counts
    current_dir = os.path.dirname(os.path.abspath(__file__))
    logging.config.fileConfig(current_dir+ "logging.conf", disable_existing_loggers=False)
    log = logging.getLogger("root")

    # gor through fastq files in the input folder
    all_fastq = glob.glob(arguments.fastq)
    for f in all_fastq:
        # read 1 is AD and read 2 is DB
        if not "_R1" in os.path.basename(f):
            # ignore R2
            continue
        ad_base = os.path.basename(f).split("_R1")[0]
        dir_name = os.path.dirname(f)
        # find DB
        db = [i for i in all_fastq if "_R2" in i and i.split("_R2")[0]==ad_base][0]
        db = os.path.join(dir_name, db)

        # process AD and DB to extract group information
        # this depends on the input mode
        AD_GROUP, DB_GROUP, AD_REF, DB_REF = parse_input_files(arguments.mode, ad_base)

        # assume all the fastq files have the filename: y/hAD*DB*_GFP_*
        output_dir_name = ad_base.split("_GFP_")[0]+"/"
        output_dir = os.path.join(arguments.o, output_dir_name)
        if not os.path.isdir(output_dir):
            os.system("mkdir -p "+output_dir)

        # make sh dir to save submission scripts
        sh_dir = os.path.join(arguments.o, "GALEN_jobs")
        if not os.path.isdir(output_dir):
            os.system("mkdir -p "+sh_dir)

        r1_csv, r2_csv = "N", "N"
        # if user wants to do alignments
        if arguments.alignment:
            r1_csv, r2_csv, sh_file = alignment.bowtie_align(f, db, AD_REF, DB_REF, output_dir, sh_dir)
            # submit job and wait
            # os.system(f"sbatch {sh_file}")
            # with r1_csv and r2_csv, add python command for read counts
            # read count script
            rc_script = os.path.join(current_dir, "read_counts.py")
            rc_cmd = f"{rc_script} -r1 {r1_csv} -r2 {r2_csv} --AD_GROUP {AD_GROUP} --DB_GROUP {DB_GROUP} --mode {arguments.mode} " \
                     f"--cutoff {arguments.cutOff} -o {output_dir}"

            with open(sh_file, "a") as f:
                f.write(rc_cmd+"\n")
            os.system(f"sbatch {sh_file}")
        else:
            # retrieve r1_csv and r2_csv
            for f in os.listdir(output_dir):
                # go through the output dir and find the csv files
                if ad_base in f:
                    if "_R1" in f and ".csv" in f:
                        r1_csv = f
                    if "_R2" in f and ".csv" in f:
                        r2_csv = f
             # check r1 and r2
            if not os.path.isfile(r1_csv) or not os.path.isfile(r2_csv):
                raise FileNotFoundError("Alignment script did not finish properly, check log")

            read_counts.RCmain(r1_csv, r2_csv, AD_GROUP, DB_GROUP, arguments.mode, output_dir, arguments.cutOff)




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

        AD_REF = param.yREF_PATH + "y_AD_wnull_" + AD_GROUP
        DB_REF = param.yREF_PATH + "y_DB_wnull_" + DB_GROUP

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
        AD_REF = param.hREF_PATH + "h_AD_" + AD_GROUP
        DB_REF = param.hREF_PATH + "h_DB_" + DB_GROUP

    elif mode == "virus":
        # human vs virus pairwise
        m = re.match(r"([hv])(AD([0-9]+)|ADNC|AD2u|ADall)([hv])(DB([0-9]+)|DBNC|DB2u|DBall)", ad_base)

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
            AD_REF = param.vREF_PATH + "v_" + AD_GROUP
        else:
            AD_REF = param.hvREF_PATH + "h_" + AD_GROUP

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

        AD_REF = param.hvREF_PATH + "h_AD_wnull_" + AD_GROUP
        DB_REF = param.heREF_PATH + "h_DB_" + DB_GROUP

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

    parser.add_argument("--alignment", action="store_true", help= "turn on alignment")
    parser.add_argument("--readCount", action="store_true", help= "turn on read counting")

    parser.add_argument("-r", help= "redirect output files to this dir")
    parser.add_argument("--cutOff", type=int, help = "assign cut off", default=20)

    args = parser.parse_args()

