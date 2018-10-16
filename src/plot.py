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


def norm_score_corr(sample_name, s, s_prime):

    score = list(s.values.flatten())
    norm_score = list(s_prime.values.flatten())
    
    pcc = pearsonr(score, norm_score)
    plt.plot(score, norm_score, ".")
    plt.xlabel("Score")
    plt.ylabel("Normalized Score")

    plt.text(0, 9, "pcc:"+str(round(pcc[0], 2)))
    plt.title(sample_name+" corr (log)", fontsize= 16)
    plt.savefig(sample_name+"_corr.png")
    plt.close()


def freq_corr(freq_one, freq_two):
    # convert matrix to list
    freq_one_list = list(freq_one.values.flatten())
    freq_two_list = list(freq_two.values.flatten())
    
    pcc = pearsonr(freq_two_list, freq_one_list)
    # plot corr
    plt.plot(np.log(freq_one_list), np.log(freq_two_list), ".")
    plt.text(-16, -5, "pcc:"+str(round(pcc[0], 2)))
    plt.savefig("test_freq_corr.png")
    plt.close()


def plot_diff(diff_list):
    x = range(len(diff_list))
    plt.plot(x, diff_list)
    plt.savefig("test_diff.png")
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
