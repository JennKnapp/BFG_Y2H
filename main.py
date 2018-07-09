import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # if no summary file provided, skip this step
    parser.add_argument('--summary', type=str, help="Summary file for making referece fasta")
    
    # path to reference file
    parser.add_argument("--ref", type=str, help="Path to reference fasta")
