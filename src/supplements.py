
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


