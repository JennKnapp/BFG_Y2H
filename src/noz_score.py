from __future__ import division
import pandas as pd
import numpy as np
import os
import datetime
from plot import *

# calculating scores based on nozomu's paper
# load the file in the same way as score.py
# detail about how to calculate the score on google doc

def pre_freq(rc_pre):
    """
    frequency for GFP pre. same as non-selective condition 
    in Nozomu's paper: use marginal frequencies
    """
    total_reads = rc_pre.values.sum()

    col_freq = rc_pre.sum(axis=0)/total_reads
    row_freq = rc_pre.sum(axis=1)/total_reads
    return row_freq, col_freq


def freq(rc_matrix):
    """
    freq = read counts for this pair / total reads in this matrix
    """
    total_reads = rc_matrix.values.sum()
    # normalize the score by +1
    rc_matrix += 1
    freq = rc_matrix / total_reads
    return freq


def get_score(pre_freq, med_freq, high_freq):
    
    """
    calculate raw score
    """
    # convert pre_freq to dataframe
    pre_freq = pd.DataFrame(data = pre_freq, columns=med_freq.columns.tolist(), index=med_freq.index.tolist())

    s = (med_freq+high_freq)/pre_freq

    return s

def norm_score(s, q):
    """
    Normalize the raw score matrix s
    q: float, quantile 
    """
    # get median of all DB scores 
    med = s.median(axis=0)
    beta = s[s>med]-med
    
    beta_q = beta.quantile(q)
    
    # if s_ij - med(s_i) < beta_i
    m1 = (s-med) < beta_q
    # if s_ij - med(s_i) >= beta_i
    m2 = (s-med) >= beta_q

    s[m1] = 1
    s[m2] = (s-med)/beta_q
    
    return s


def get_rank(norm_s):
    """
    In total we have 4 ranks for each protein pair
    aBC1-bBC1, aBC1-bBC2, aBC2-bBC1, aBC2-bBC2 
    """
    print datetime.datetime.now()
    # unstack the matrix
    transform = norm_s.unstack().reset_index()
    # rename col
    transform.columns = ['DB', 'AD', 's_prime']
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB'].str.split('_', expand=True)
    transform[['AD', 'AD_BC']] = transform['AD'].str.split('_', expand=True)
    
    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    print len(transform.groupby(['Interaction']).groups.keys())
    print datetime.datetime.now()
    return transform

def main(GFP_pre, GFP_med, GFP_high):

    #calculate GFP_pre freq
    row_freq, col_freq = pre_freq(GFP_pre)
    GFP_pre_freq = np.outer(row_freq, col_freq)

    med_freq = freq(GFP_med)
    high_freq = freq(GFP_high)

    # test plot
    # should be very little corr
    # freq_corr(high_freq, med_freq)
    
    # get raw scores
    s = get_score(GFP_pre_freq, med_freq, high_freq)
    test = s.values.flatten()
    test.sort()
    plot_diff(test)
    
    percentile = np.arange(0.1, 1.0, 0.05)
    norm_s = norm_score(s, 0.75)
    rank = get_rank(norm_s)

    rank = rank.sort_values(by=['s_prime'], ascending=False)
    print rank
    plot_s(rank, "s_prime_sorted.png")

#    print norm_s
    # test different rho
#    for p in percentile:
        # get normalized scores
#        norm_s = norm_score(s, p)
#        sample_name = "rho = "+ str(p)
        # test corr of score and normed score
#        norm_score_corr(sample_name, s, norm_s)



if __name__ == "__main__":
    # test on yAD2 DB1
    test_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/rerun_analysis/yAD4DB1/"
    
    os.chdir(test_dir)
    for f in os.listdir(test_dir):
        if not f.endswith("_combined_counts.csv"):
            continue
        fname = "./"+f
        if "pre" in f:
            GFP_pre = pd.read_table(fname, sep=",", index_col=0)
        elif "med" in f:
            GFP_med = pd.read_table(fname, sep =",", index_col=0)
        elif "high" in f:
            GFP_high = pd.read_table(fname, sep =",", index_col=0)

    main(GFP_pre, GFP_med, GFP_high)            


