import pandas as pd
import os


def read_mcc(output_dir):
    """
    read mcc summary from directory and make a summary for each 
    group
    output format:
    group   weight  floor   rank    max_mcc
    """

    frames= []
    for d in os.listdir(output_dir):
        for f in os.listdir(output_dir+d):
            if f == "max_parameters.csv":
                p = pd.read_csv(output_dir+d+"/"+f)
                p["group"] = d
                
                # for dk
                dk_mcc_yi1 = output_dir+d+"/DK_mcc_summary_yi1.csv"
                df = pd.read_csv(dk_mcc_yi1)
                # select rows with max mcc
                max_mcc_summary = df[(df["floor"] == p["max_floor"].tolist()[0]) & (df["rank"] == p["max_rank"].tolist()[0]) & (df["weight"] == p["max_weight"].tolist()[0])]
                p["dk_max_mcc"] = max_mcc_summary["mcc"].max()
                #,max_rank,max_weight,max_rho,noz_max_rank
                
                # for noz
                noz_mcc_yi1 = output_dir +d+"/noz_mcc_summary_yi1.csv"
                df = pd.read_csv(noz_mcc_yi1)
                max_mcc_summary = df[(df["rho"] == p["max_rho"].tolist()[0]) & (df["rank"] == p["noz_max_rank"].tolist()[0])]
                p["noz_max_mcc"] = max_mcc_summary["mcc"].max()
                frames.append(p)            
    groups = pd.concat(frames).reset_index(drop=True)
    groups.to_csv("all_group_mcc_summary_yi1.csv", index=False)


def read_noz_score(output_dir, gene_list):
    frames = []
    genes = pd.read_table(gene_list)

    for d in os.listdir(output_dir):
        group = d
        for index, row in genes.iterrows():
            if row['group'] == group:
                print group
                raw_f =pd.read_csv(output_dir+d+"/noz_raw_score.csv")
                norm_f =pd.read_csv(output_dir+d+"/noz_norm_score.csv")
                raw_f.columns = ['DB', 'AD', 'raw_score']
                norm_f.columns = ['DB', 'AD', 'norm_score']
                raw_scores = raw_f[(raw_f.DB.str.contains(row['DB'])) & (raw_f.AD.str.contains(row['AD']))]
                norm_scores = norm_f[(norm_f.DB.str.contains(row['DB'])) & (norm_f.AD.str.contains(row['AD']))]
                raw_scores["norm_score"] = norm_scores["norm_score"].tolist()
                frames.append(raw_scores)
        print frames
    groups = pd.concat(frames).reset_index(drop=True)
    groups.to_csv("dayag_gene_scores.csv", index=False)
        
if __name__ == "__main__":
    output_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/181109_test/"
    #read_mcc(output_dir)
    gene_list = "/home/rothlab/rli/02_dev/08_bfg_y2h/gene_list.txt"
    gene_list_v1 = "/home/rothlab/rli/02_dev/08_bfg_y2h/181119_gene_list.txt"
    read_noz_score(output_dir, gene_list_v1)
