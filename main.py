import argparse
from alignment import *
from param import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # if no summary file provided, skip this step
    parser.add_argument('--create', help="Summary file for making referece fasta", nargs=3)
    
    args = parser.parse_args()

    summary = args.create

    if summary is not None:
        for f in summary:
            if "AD" in f:
                AD_summary = f
            elif "DB" in f:
                DB_summary = f
            else:
                output = f

        create_fasta(AD_summary, DB_summary, output, group_spec=True, AD="G1", DB="G4")
