###################################

# if you want to create fasta file from summary
MAKE_FASTA = False

# if you want to build index
BUILD = False

# If you want to align the reads 
# Set this to True
ALIGN = False

# If you want to do the analysis
# Set this to True
ANALYSIS = True

###################################

# AD and DB group
# if it's set to G0, means all
AD_GROUP="G1"
DB_GROUP="G4"

# default size for AD and DB
# the matrix is built with col=DB_SIZE and row=AD_SIZE
# you can change this in read_counts.

AD_SIZE = 2491
DB_SIZE = 2510

# file path
ad_fastq = ""
db_fastq = ""

# ref path

AD_REF = "/home/rothlab/rli/02_dev/08_bfg_y2h/ref/y_AD_G1"
DB_REF = "/home/rothlab/rli/02_dev/08_bfg_y2h/ref/y_DB_G4"


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
