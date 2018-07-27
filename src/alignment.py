from param import *
import os

def bowtie_align(fastq, ref, output):
    """
    Align r1 and r2 to reference
    Log bowtie output 
    """

    basename = os.path.basename(fastq)

    if "R1" in fastq:
        param = "-q --no-sq --norc --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_AD_BC.sam')
    elif "R2" in fastq:
        param = "-q --no-sq --nofw --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_DB_BC.sam')

    input_f = "-x " + ref + " -U " + fastq + " -S " + os.path.join(output, sam_file)
    
    log_f = os.path.join(output, sam_file.replace(".sam", "_bowtie.log"))
    
    command = BOWTIE2 + param + input_f + " 2> " + log_f

    os.system(command)
