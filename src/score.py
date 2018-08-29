from __future__ import division
import os
import math
import numpy as np
import pandas as pd
import logging
import logging.config
import itertools

#scores_log = logging.getLogger("scores")

def freq(matrix):
        
    # sum of matrix
    total = matrix.values.sum()

    freq_df = matrix / total

    return freq_df


def marginal_freq(matrix):
    # sum of matrix
    total = matrix.values.sum()

    row_freq = matrix.sum(axis=0)/total
    col_freq = matrix.sum(axis=1)/total

    return row_freq, col_freq


def load_YI1(yi):
    "Load gold standard for yi1"
    yi = pd.read_table(yi, sep="\t")
    AD_GOLD = yi.AD_ORF_ID
    DB_GOLD = yi.DB_ORF_ID

    #for ad in yi.AD_ORF_ID:
    #    AD_GOLD.append(ad+"_BC-1")
    #    AD_GOLD.append(ad+"_BC-2")

    #for db in yi.DB_ORF_ID:
    #    DB_GOLD.append(db+"_BC-1")
    #    DB_GOLD.append(db+"_BC-2")

    return AD_GOLD, DB_GOLD

def rename(df):
    """remove _BC-* from orf names"""
    return df.rename(columns = lambda x : str(x)[:-5], index = lambda x : str(x)[:-5])


def mcc(l, range_test):
    PRMCC = []
    # each list consist of sorted scores
    for i in range(range_test):
        temp_screen = l[i:]
        temp_ns = l[i+1:]

        precision = sum(temp_screen) / i * 100
        recall = sum(temp_screen)/sum(l) *100

        TP=sum(temp_screen)

        FN=sum(temp_ns)


        MCC= (TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))*100;
        PRMCC.append([precision, recall, MCC])

    return PRMCC


def score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, opt_mix, yi_GOLD):
    # to find optimum set of parameters
    # we test all possible combinations
    # generates a list of tuples containing all possible combinations
    l = [weights, floor_perc, opt_mix]
    comb = list(itertools.product(*l))
    
    # for GFP_pre
    row_freq, col_freq = marginal_freq(GFP_pre)
    
    # for GFP_med
    med_freq = rename(freq(GFP_med)) 
    
    # for GFP_high
    high_freq = rename(freq(GFP_high))

    shape = GFP_pre.shape
    total_rows = shape[0] #AD
    total_cols = shape[1] #DB

    row_sorted = sorted(row_freq.tolist())
    col_sorted = sorted(col_freq.tolist())

    # get gold standard
    AD_GOLD, DB_GOLD = load_YI1(yi_gold)

    # get ad and db genes
    
    AD_NAMES = list(med_freq.index) # rownames
    DB_NAMES = list(med_freq) #colnames

    # find intersection
    AD_intersect = AD_GOLD and AD_NAMES
    DB_intersect = DB_GOLD and DB_NAMES

    for s in comb:
        weight = s[0]
        floor = s[1]
        mix_index = int(s[2])

        # log
        #scores_log.info("iteration: w - %d, f - %d, m - %d", weight, floor, mix)
        
        row_cut = row_sorted[round(total_rows/floor)]
        col_cut = col_sorted[round(total_cols/floor)]

        AD_freq = row_freq.where(row_freq > row_cut, row_cut)
        DB_freq = col_freq.where(col_freq > col_cut, col_cut)

        # rebuild matrix from two vectors
        freq_mx = np.outer(AD_freq, DB_freq)

        # convert to df with index and row names 
        pre_freq = rename(pd.DataFrame(data = freq_mx, columns = DB_freq.index.tolist(), index = AD_freq.index.tolist()))
        
        # log
        #scores_log.info("Shape of GFP_pre ")
        # score
        IS = ((weight * high_freq) + med_freq) / pre_freq
        
        # plot scores (debug)

        # replace with nan
        IS = IS.replace(0, np.nan)
        
        # log2 median of DB
        log_med = np.log2(IS.median(axis=0))
        
        # log2 of matrix
        log_is = np.log2(IS)

        # subtract
        IS_normed = log_is.sub(log_med)
        
        # plot normalized scores (debug)
        
        # test which set of bc we want to use
        # in this case, we will have max 4 min 1 score(s) for each orf pair
        interaction_score = []
        
        all_pairs = list(itertools.product(AD_intersect, DB_intersect))
        for pair in all_pairs:
            # 0 is AD ; 1 is DB
            # get this pair from IS_normed
            scores = IS_normed.loc[pair[0], pair[1]]
            fla = [val for sublist in scores.values.tolist() for val in sublist if str(val) !='nan']
            fla.sort()
            try:
                s = fla[mix_index]
            except Exception e:
                s = np.nan
            if str(s) == 'nan': 
                continue
            else:
                interaction_score.append([pair[0], pair[1], s])

        sort_is = pd.DataFrame(interaction_score, columns=["AD", "DB", "Score"])
        sort_is[["Score"]] = sort_is[["Score"]].apply(pd.to_numeric)
        
        
        sort_is = sort_is.sort_values(by="Score", ascending=False)
        
        # select intersection of two dfs


if __name__ == "__main__":
    weights = np.arange(0, 2, 0.2)
    floor_perc = np.arange(5,20,1)
    opt_mix = [1,2,3,4]
    
    GFP_pre = "" 
    GFP_med = ""
    GFP_high = ""
    
    score_main(weights, floor_perc,opt_mix)
    


