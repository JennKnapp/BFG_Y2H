###################################

# if you want to create fasta file from summary
MAKE_FASTA = False

# If you want to do alignments
# Set this to True
ALIGNMENT = True

# if you want to do read counts
READ_COUNT = False

###################################

#AD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_AD.csv"
#DB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_DB_AA.csv"

AD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_AD.csv"
DB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_DB_AA.csv"


###################################

# program path

BOWTIE2 = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2 "
BOWTIE2_BUILD = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build "
SAMTOOLS = "/home/rothlab/rli/lib/samtools-1.4.1/samtools "

###################################

# Padding sequences 
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
