###################################

# If you want to align the reads 
# Set this to True
ALIGN = True

# If you want to do the analysis
# Set this to True
ANALYSIS = False

###################################

# default size for AD and DB
# the matrix is built with col=DB_SIZE and row=AD_SIZE
# you can change this in read_counts.

AD_SIZE = 2491
DB_SIZE = 2510

# file path

ad_fastq = ""
db_fastq = ""

# ref path

AD_REF = "/home/rothlab/rli/02_dev/08_bfg_y2h/ds_ref/ds_AD_ref"
DB_REF = "/home/rothlab/rli/02_dev/08_bfg_y2h/ds_ref/ds_DB_ref"


###################################

# program path

BOWTIE2 = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2 "
BOWTIE2_BUILD = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build "

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
