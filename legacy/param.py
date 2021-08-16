#!/usr/bin/env python3.7
import numpy as np

###################################

# If you want to do alignments
# Set this to True
ALIGNMENT = True

# if you want to do read counts
# Set this to True
READ_COUNT = True
cut_off = 20 

##################################

# for score optimization

# gold standard to use
yGOLD="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/YI_1.txt"
hGOLD="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/FIPP.txt"
hGOLD="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/FIPP_0.6.txt"
# parameters to test DK's normalization method
#weights = np.arange(0, 2.4, 0.2)

pho = np.arange(0.1, 0.8, 0.05)
floor_perc = np.arange(5,25,2.5)
weights = [1.0]
mix_index= [0,1,2,3]

#floor_perc = [10.0]
#floor_perc = np.arange(5,10,5)
#weights = [1.0]
#mix_index= [0,1]
#pho = np.arange(0.1, 0.2, 0.05)

# Our ORFs are double barcoded and 
# in this case we test all 4 indexes

# for calculating MCC
hlitBM = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/FIPP.txt"
hlitBM = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/FIPP_0.6.txt"
ylitBM13="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/litbm_13.txt"

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

# program path 
# change according to your system

# BOWTIE2 = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2 "
# BOWTIE2_BUILD = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build "
# SAMTOOLS = "/home/rothlab/rli/lib/samtools-1.4.1/samtools "

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
