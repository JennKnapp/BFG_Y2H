import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr

def bc_corr(sample_name, uptag_matrix, dntag_matrix):
    # convert matrix to list
    up_tag_list = list(uptag_matrix.values.flatten())
    dn_tag_list = list(dntag_matrix.values.flatten())
    pcc = pearsonr(up_tag_list, dn_tag_list)
    # plot corr
    # x as up tags
    # y as dn tags
    plt.plot(np.log(up_tag_list), np.log(dn_tag_list), ".")
    plt.xlabel('uptag', fontsize=14)
    plt.ylabel('dntag', fontsize=14)
    plt.ylim([0, 15])
    plt.xlim([0,15])

    plt.text(1, 14, "pcc:"+str(round(pcc[0], 2))+" RC:"+str(sum(up_tag_list)))

    plt.title(sample_name+"_bc_corr (log)", fontsize= 16)
    plt.savefig(sample_name+"_bc_corr.png")
    plt.close()

if __name__ == "__main__":
    sample_name = "yAD1DB4_high" 
    uptag_matrix = "/home/rothlab/rli/02_dev/08_bfg_y2h/allbyall_output/yAD1DB4_high/uptag_rawcounts.csv"
    dntag_matrix = "/home/rothlab/rli/02_dev/08_bfg_y2h/allbyall_output/yAD1DB4_high/dntag_rawcounts.csv"

    uptag_matrix = pd.read_csv(uptag_matrix, index_col=0).astype(int)
    dntag_matrix = pd.read_csv(dntag_matrix, index_col=0).astype(int)

    bc_corr(sample_name, uptag_matrix, dntag_matrix)
