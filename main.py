import argparse
from create_fasta import *
from param import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # if no summary file provided, skip this step
    parser.add_argument('--create', help="Summary file for making referece fasta", nargs=3)
    parser.add_argument('--fasta'. help="Path to fasta file")

    args = parser.parse_args()
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
        fasta_output = args.fasta

    if fasta_output is not None:
        list_fasta = os.listdir(fasta_output)
        for fasta in list_fasta:
            build_index(os.path.join(output,fasta), output)

    
