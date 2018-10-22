from __future__ import division
import os
import argparse
import datetime
import math
import numpy as np
import pandas as pd
import logging
import logging.config
import itertools
from plot import *
from param import *
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef

def freq(matrix):
        
    # sum of matrix
    total = matrix.values.sum()
    freq_df = matrix / total

    return freq_df


def marginal_freq(matrix):
    # sum of matrix
    total = matrix.values.sum()
    col_freq = matrix.sum(axis=0)/total
    row_freq = matrix.sum(axis=1)/total
    return row_freq, col_freq


def load_YI1(yi):
    "Load gold standard for yi1"
    yi = pd.read_table(yi, sep="\t")
    yi.columns = ["AD", "DB"]
    yi["Interactions"] = yi[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return yi


def rename(df):
    """remove _BC-* from orf names"""
    return df.rename(columns = lambda x : str(x)[:-5], index = lambda x : str(x)[:-5])


def calculate_freq(GFP_pre, GFP_high, GFP_med):
    
    # for GFP_pre
    row_freq, col_freq = marginal_freq(GFP_pre)
    # for GFP_med
    med_freq = freq(GFP_med)
    # for GFP_high
    high_freq = freq(GFP_high)

    med_freq_renamed = rename(freq(GFP_med))

    AD_NAMES = list(med_freq_renamed.index)
    DB_NAMES = list(med_freq_renamed)

    #heat_freq(med_freq)
    
    return row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES


def get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut):
    # dicts contains intereactions and scores for each index 
    # compare to gold standard and plot PRC
    mcc = {}
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()
    for j in dicts.keys():
        i = dicts[j]
        #i["Interaction"] = i.index
        i = i.reset_index()
        # compare gold with our set
        i["ishit_interaction"] = i.Interaction.isin(gold_st.Interactions).astype(int)
        
        ## build screen set
        # AD has to appear in AD_GOLD at least once
        # DB has to appear in DB_GOLD at least once
        i["is_AD"] = i.AD_x.isin(AD_GOLD)
        i["is_DB"] = i.DB_x.isin(DB_GOLD)
         
        # frequencies should be greater than floow
        i = pd.merge(i, AD_freq, on="AD_name")
        i = pd.merge(i, DB_freq, on="DB_name")
        
        i.rename(columns={"0_x":"AD_mfreq", "0_y":"DB_mfreq"}, inplace=True)
        i["AD_cut"] = i.AD_mfreq > row_cut
        i["DB_cut"] = i.DB_mfreq > col_cut
        # get union
        
        i["screen"] = (i.is_AD & i.is_DB & i.AD_cut & i.DB_cut).astype(int)
        i = i.sort_values(by="Score", ascending=False)  
        # hit
        y_pred = i.ishit_interaction.tolist()
        # screen 
        y_network = i.loc[i.screen == 1].ishit_interaction.tolist()
        # scores
        y_score = i.Score.tolist()
        MAXMCC = prcmcc(y_network, 1000)
        mcc[j] = MAXMCC
    return mcc


def prcmcc(label, test_range):
    # pred: predicted labels
    # label: actual labels
    PRCMCC = []
    total_screen = len(label)
    if total_screen < test_range:
        test_range = total_screen-1
    
    for i in range(1,test_range+1):

        test_screen = label[:i]
        test_nonscreen = label[i:]
        # true positive: condition positive and predicted pos
        TP = sum(test_screen)
        # true negative: condition negative and predicted neg
        TN = test_nonscreen.count(0)
        # false positive: condition negative and predicted pos
        FP = test_screen.count(0)
        # false negative: condition positive and predicted neg
        FN = sum(test_nonscreen)
        # precision = TP / (TP+FP)
        precision = sum(test_screen)/i * 100
        # recall = TP / (TP+FN)
        recall = sum(test_screen)/sum(label) * 100
        MCC= (TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))*100 
        PRCMCC.append([precision, recall, MCC])
    
    df = pd.DataFrame(PRCMCC, columns=["precision", "recall", "MCC"])

    MAXid = df.MCC.idxmax()
    #print MAXid
    MAXMCC = df.loc[MAXid].tolist()
    #all_mcc.append(MAXMCC) 
    MAXMCC.insert(0, total_screen)
    return MAXMCC


def test_rank(IS_normed, all_pairs, mix_index):

    dicts = {"is_{}".format(i): {} for i in mix_index}
    all_pairs = pd.DataFrame(all_pairs, columns=['AD', 'DB'])
    all_pairs['Interaction'] = all_pairs.AD.str.cat(all_pairs.DB, sep="_")

    transform = IS_normed.unstack().reset_index()
    transform.columns = ['DB_name', 'AD_name', 'Score']
    # drop nan score
    #transform = transform.dropna(how='any')
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB_name'].str.split('_', expand=True)
    transform[['AD', 'AD_BC']] = transform['AD_name'].str.split('_', expand=True)
    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    merged = pd.merge(all_pairs, transform, on="Interaction")
    g = merged.sort_values(["Score"], ascending=False).groupby(['Interaction'])
    

    for i in mix_index:
        d_name = "is_{}".format(i)
        dicts[d_name] = g.nth(i).dropna(how='any')
    
    return dicts


def score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, gold_st):
        
    # to find optimum set of parameters
    # we test all possible combinations
    # generates a list of tuples containing all possible combinations
    l = [weights, floor_perc]
    comb = list(itertools.product(*l))
    
    # get frequencies

    row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES = calculate_freq(GFP_pre, GFP_high, GFP_med)
    
    shape = GFP_pre.shape
    total_rows = shape[0] #AD
    total_cols = shape[1] #DB

    row_sorted = sorted(row_freq.tolist())
    col_sorted = sorted(col_freq.tolist())
    
    # get gold standard
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()

    # find intersection
    AD_intersect = list(set(AD_GOLD) & set(AD_NAMES))
    DB_intersect = list(set(DB_GOLD) & set(DB_NAMES))
    
    all_pairs = list(itertools.product(AD_intersect, DB_intersect))
    
    mix_index = [0,1,2,3]
    
    # optimization
 
    output_csv = {}
    for s in comb:
        weight = s[0]
        floor = s[1]
        
        row_cut = row_sorted[int(round(total_rows/floor))]
        col_cut = col_sorted[int(round(total_cols/floor))]

        AD_freq = row_freq.where(row_freq > row_cut, row_cut)
        DB_freq = col_freq.where(col_freq > col_cut, col_cut)

        # rebuild matrix from two vectors
        freq_mx = np.outer(AD_freq, DB_freq)
        
        pre_freq = pd.DataFrame(data = freq_mx, columns = DB_freq.index.tolist(), index = AD_freq.index.tolist())
    
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
        
        dicts = test_rank(IS_normed, all_pairs, mix_index)
        
        AD_freq = AD_freq.to_frame()
        DB_freq = DB_freq.to_frame()
        
        AD_freq["AD_name"] = AD_freq.index
        DB_freq["DB_name"] = DB_freq.index
        
        output_name = ""
        mcc_list = get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut)
        output_csv[(weight, floor)] = mcc_list
    return output_csv

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='BFG-Y2H Scores')
    parser.add_argument("--sample", help="Path for read counts files")

    args = parser.parse_args()
    counts = args.sample

    sample_name = os.path.basename(counts)    
    
    os.chdir(counts)
    # find the optimized parameters
    yi = load_YI1(GOLD)
    for f in os.listdir(counts):
        if not f.endswith("_combined_counts.csv"):
            continue
        fname = "./"+f        
        if "pre" in f:
            GFP_pre = pd.read_table(fname, sep=",", index_col=0)
        elif "med" in f:
            GFP_med = pd.read_table(fname, sep =",", index_col=0)
        elif "high" in f:
            GFP_high = pd.read_table(fname, sep =",", index_col=0)
    
    output_csv = score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, yi)
    df = pd.DataFrame(output_csv).unstack().reset_index()
    df.columns = ["weight", "floor","index", "mcc_list"]
    df[["total_screen", "precision","recall","mcc"]] = pd.DataFrame(df.mcc_list.tolist(), columns=["total_screen", "precision","recall","mcc"])
    df = df.drop(columns=["mcc_list"])
    df["sample"] = sample_name
    df.to_csv("./score_optimization_all.csv", index=False)
    
    # get the top 3 MCC for this sample and save to a different file

    maxdf = df.nlargest(5, "mcc")
    maxdf.to_csv("./score_optimization_max.csv", index=False)
