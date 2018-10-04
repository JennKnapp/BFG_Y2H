import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
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

    plt.title(sample_name+"bc_corr (log)", fontsize= 16)
    plt.savefig(sample_name+"bc_corr.png")
    plt.close()

def bar_freq(pre_freq, med_freq, high_freq):
    # convert matrix to list and remove zeros
    print "Making bar plot"
    all_freq = [pre_freq, med_freq, high_freq]
    
    for i in all_freq:
        l = i.values.tolist()
        l = reduce(lambda x,y:x+y,l)
        print len(l)

        l.sort()
        l = l[:1000]
        
        # make a bar plot
        plt.bar(np.arange(len(l)), l, align="center")
        plt.title("Freq")
        plt.savefig("/home/rothlab/rli/www/html/bar_freq.png")
        plt.close()
        print "done"
        break

def heat_freq(freq_matrix):
    
    hm = sns.heatmap(freq_matrix, xticklabels= False, yticklabels=False, annot=False)
    hm = hm.get_figure()
    hm.savefig("/home/rothlab/rli/www/html/heatmap_output.png")
    

def plot_prc(precision, recall, output_file):

    plt.step(recall, precision, color='b', alpha=0.6,where='post')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.savefig(output_file)
    plt.close()


def rank_prc(precision, recall, output_file):
    #x = len(precision) + 10
    #y = range(0, 100)
    plt.plot(precision, color= "b", alpha=0.6)
    plt.plot(recall, color="r", alpha=0.6)
    plt.xlabel("rank")
    plt.ylabel("percentage")
    plt.ylim((0,100))
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    sample_name = "yAD1DB4_high" 
    uptag_matrix = "/home/rothlab/rli/02_dev/08_bfg_y2h/allbyall_output/yAD1DB4_high/uptag_rawcounts.csv"
    dntag_matrix = "/home/rothlab/rli/02_dev/08_bfg_y2h/allbyall_output/yAD1DB4_high/dntag_rawcounts.csv"

    uptag_matrix = pd.read_csv(uptag_matrix, index_col=0).astype(int)
    dntag_matrix = pd.read_csv(dntag_matrix, index_col=0).astype(int)

    bc_corr(sample_name, uptag_matrix, dntag_matrix)
