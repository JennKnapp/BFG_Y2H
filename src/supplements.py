import pandas as pd

def read_summary(AD_sum, DB_sum, AD_group="G0", DB_group="G0"):
    """
    Read from AD and DB summary files. 
    Grep gene names based on AD and DB group
    If group == G0, grep all
    """
    
    AD_summary = pd.read_table(AD_sum, sep="\t")
    DB_summary = pd.read_table(DB_sum, sep="\t")

    # grep group 
    if AD_group!="G0":
        AD_summary = AD_summary[AD_summary.Group==AD_group]
    if DB_group!="G0":
        DB_summary = DB_summary[DB_summary.Group==DB_group]
    
    # grep gene names
    AD_genes = AD_summary.Locus.tolist()
    DB_genes = DB_summary.Locus.tolist()

    return AD_genes, DB_genes


def parse_ds_ref(fasta):
    """
    separate dayag's fasta file
    """

    with open(fasta, "r") as ref, open("ds_AD_ref.fasta", "w") as ad, open("ds_DB_ref.fasta", "w") as db:
        c = ref.readlines()
        line =0
        while line in range(len(c)):
        #for line in c:
            print line
            if ">c" in c[line]: 
                line+=2 
                continue
            
            if ">AD" in c[line]:
                ad.write(c[line])
                ad.write(c[line+1])
            else:
                db.write(c[line])
                db.write(c[line+1])
            line+=2

if __name__ == "__main__":
    fasta = "./ds_ref/barcodes.fasta"
    parse_ds_ref(fasta)


