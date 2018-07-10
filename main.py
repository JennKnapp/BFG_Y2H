import argparse
from create_fasta import *
from param import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # if no summary file provided, skip this step
    parser.add_argument('--create', help="Summary file for making referece fasta", nargs=3)
    parser.add_argument('--buildindex'. help="Path to fasta file")
    
    parser.add_argument("--ad", help="Read one and corresponding reference file", nargs=2)
    parser.add_argument("--db", help="Read two and corresponding reference file", nargs=2)
    
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
        
        create_fasta(AD_summary, DB_summary, fasta_output, group_spec=True, AD="G1", DB="G4")

    if fasta_output is None:
        fasta_output = args.buildindex

    if fasta_output is not None:
        list_fasta = os.listdir(fasta_output)
        for fasta in list_fasta:
            build_index(os.path.join(output,fasta), output)

    
    # Alignment 
    
    if ALIGN:
        ad = args.ad
        db = args.db
        sam_output = args.output

        bowtie_align(ad[0], ad[1], sam_output)
        bowtie_align(db[0], db[1], sam_output)

    # Read counts

