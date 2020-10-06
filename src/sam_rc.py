import argparse
import pandas as pd
import numpy as np
import param
import os 
import supplements

def preprocess_sam(r1_sam, r2_sam):
    """
    preprocess sam files 
    
    """
    #if not os.path.isfile(r1_sam) or not os.path.isfile(r2_sam):
    #    print("file doesn't exist")
    #    exit(0)

    dir_name = os.path.dirname(r1_sam)
    r1_basename = os.path.basename(r1_sam)
    r2_basename = os.path.basename(r2_sam)

    sorted_r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_sorted.sam"))
    
    sort_r1 = param.SAMTOOLS+"sort -n -o "+sorted_r1+" "+r1_sam
    sorted_r2 = os.path.join(dir_name, r2_basename.replace(".sam","_sorted.sam"))
    sort_r2 = param.SAMTOOLS+"sort -n -o "+sorted_r2+" "+r2_sam

    # remove headers
    r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_noh.sam"))
    r2 = os.path.join(dir_name, r2_basename.replace(".sam", "_noh.sam"))
    os.system(sort_r1)
    os.system(sort_r2)

    #os.system("rm "+r1_sam)
    #os.system("rm "+r2_sam)

    os.system("grep -v \"^@\" "+sorted_r1+" > "+r1)
    os.system("grep -v \"^@\" "+sorted_r2+" > "+r2)
    r1_csv = os.path.join(dir_name, r1.replace(".sam", ".csv"))
    r2_csv = os.path.join(dir_name, r2.replace(".sam", ".csv"))
    os.system("cut -f 1-5 "+r1+" > "+ r1_csv)
    os.system("cut -f 1-5 "+r2+" > "+ r2_csv)
    os.system("rm "+r1)
    os.system("rm "+r2)

    return r1_csv, r2_csv
    
def read_count_hap(r1_csv, r2_csv, DB_genes):
    empty_matrix = pd.DataFrame(0, index = DB_genes, columns = DB_genes)
    
    f1 = open(r1_csv, "rb")
    f2 = open(r2_csv, "rb")

    i = True
    lines = 0
    pairs = {}

    fail = 0
    count = 0
    while 1:
        r1_line = f1.readline() # uptag
        r2_line = f2.readline() # dntag

        if r1_line == "" or r2_line == "":
            i = False
            print("End of file")
            break

        r1_line = r1_line.strip().split("\t")
        r2_line = r2_line.strip().split("\t")
        
        if r1_line[0] != r2_line[0]:
            i = False
            print("# READ ID DOES NOT MATCH #")
            break
        
        if int(r1_line[4]) < param.cut_off or int(r2_line[4]) < param.cut_off: # check quality
            fail += 1
            continue

        if r1_line[2] == "*" or r2_line[2] =="*":
            fail +=1
            continue

        r1_name = r1_line[2].split(";")
        r2_name = r2_line[2].split(";")

        if r1_name[-1] != r2_name[-1]:
            count+=1
            pairs[(r2_name[1], r1_name[1])] = pairs.get((r2_name[1], r1_name[1]), 0) + 1
    matrix = (pd.Series(pairs)
            .unstack(fill_value=0)
            .T
            .reindex(index=empty_matrix.index, columns=empty_matrix.columns, fill_value=0))
    f1.close()
    f2.close()
    
    diag = pd.Series(np.diag(matrix), index=[matrix.index, matrix.columns])
    print(diag)
    
    return diag

def read_DB(hDB):
    """
    get a list of db gene from hDB summary
    """
    summary = pd.read_table(hDB, sep="\t")

    DB_genes = summary.Locus.tolist()
    return DB_genes

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='read count from sam files')
    
    parser.add_argument("--r1sam", help="sam file for read one")
    parser.add_argument("--r2sam", help="sam file for read two")
    parser.add_argument("--mode", help="human or yeast")

    parser.add_argument("--r1csv", help="csv file for read one")
    parser.add_argument("--r2csv", help="csv file for read two")
    
    parser.add_argument("--prefix", help= "output prefix")

    args = parser.parse_args()

    prefix = args.prefix
    r1_sam = args.r1sam
    r2_sam = args.r2sam

    if r1_sam:
        r1_csv, r2_csv = preprocess_sam(r1_sam, r2_sam)
        DB_genes = read_DB(param.hDB_summary) 
        diag = read_count_hap(r1_csv, r2_csv, DB_genes)
    else:

        r1_csv = args.r1csv
        r2_csv = args.r2csv
        DB_genes = read_DB(param.hDB_summary)
        diag = read_count_hap(r1_csv, r2_csv, DB_genes)
    
    diag.to_csv(prefix+"_matrix.csv")
