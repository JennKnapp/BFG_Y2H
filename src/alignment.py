import param
import os
import logging
import logging.config
import argparse

#logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
align_log = logging.getLogger("alignments")

def bowtie_align(fastq, ref, output):
    """
    Align r1 and r2 to reference
    Log bowtie output 
    """
    align_log.info("Starts aligning... %s", fastq)
    align_log.info("Reference: %s", ref)
    align_log.info("Output: %s", output)

    basename = os.path.basename(fastq)

    if "R1" in fastq:
        params = "-q --norc --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_AD_BC.sam')
    elif "R2" in fastq:
        
        params = "-q --nofw --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_DB_BC.sam')

    input_f = "-x " + ref + " -U " + fastq + " -S " + os.path.join(output, sam_file)
    
    log_f = os.path.join(output, sam_file.replace(".sam", "_bowtie.log"))
    
    command = param.BOWTIE2 + params + input_f + " 2> " + log_f

    os.system(command)
    print(command)

    align_log.info("Alignment finished for %s", fastq)

    return os.path.join(output, sam_file)


def bowtie_align_hap(fastq, ref, output):

    """
    align hDB to all the hDB ref
    align R1 to hDB uptag and R2 to hDB downtag
    """

    basename = os.path.basename(fastq)

    if "R1" in fastq:
        params = "-q --norc --local --very-sensitive-local -t -p 23 -k 2 --reorder "
        sam_file = basename.replace('.fastq.gz','_DB_BC_up.sam')
    elif "R2" in fastq:  
        params = "-q --nofw --local --very-sensitive-local -t -p 23 -k 2 --reorder "
        sam_file = basename.replace('.fastq.gz','_DB_BC_dn.sam')

    input_f = "-x " + ref + " -U " + fastq + " -S " + os.path.join(output, sam_file)
    log_f = os.path.join(output, sam_file.replace(".sam", "_bowtie.log"))
    command = param.BOWTIE2 + params + input_f + " 2> " + log_f
    os.system(command)

    return os.path.join(output, sam_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")

    args=parser.parse_args()

    f = args.fastq
    ref = "/home/rothlab/rli/02_dev/08_bfg_y2h/h_ref/h_DB_all"
    output = "/home/rothlab/rli/02_dev/08_bfg_y2h/output/190821_hDB/hDB_test_bowtie_k/"

    for fastq in os.listdir(f):
        fastq = f+fastq
        bowtie_align_hap(fastq, ref, output)

