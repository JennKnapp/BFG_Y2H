import os
import argparse
import datetime
import math
import numpy as np
import pandas as pd
import logging
import logging.config
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef
#scores_log = logging.getLogger("scores")

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
    AD_GOLD = yi.AD_ORF_ID.tolist()
    DB_GOLD = yi.DB_ORF_ID.tolist()

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


def plot_prc(precision, recall, output_file):

    plt.step(recall, precision, color='b', alpha=0.6,where='post')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.savefig(output_file)
    plt.close()





def score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, yi_gold, output):
    
    # output directory
    d_name = os.path.join(output, 
    
    
    
    
    # to find optimum set of parameters
    # we test all possible combinations
    # generates a list of tuples containing all possible combinations
    l = [weights, floor_perc]
    comb = list(itertools.product(*l))
    
    # for GFP_pre
    row_freq, col_freq = marginal_freq(GFP_pre)
   
    # for GFP_med
    med_freq_renamed = rename(freq(GFP_med)) 
    med_freq = freq(GFP_med)
    # for GFP_high
    #high_freq = rename(freq(GFP_high))
    high_freq = freq(GFP_high)
    
    shape = GFP_pre.shape
    total_rows = shape[0] #AD
    total_cols = shape[1] #DB

    row_sorted = sorted(row_freq.tolist())
    col_sorted = sorted(col_freq.tolist())
    
    # get gold standard
    AD_GOLD, DB_GOLD = load_YI1(yi_gold)
    # create df for gold standard
    gold_st = pd.DataFrame({"AD":AD_GOLD, "DB":DB_GOLD})
    gold_st["Interactions"] = gold_st[['AD', 'DB']].apply(lambda x: '-'.join(x), axis=1)
    # get ad and db genes
    
    AD_NAMES = list(med_freq_renamed.index) # rownames
    DB_NAMES = list(med_freq_renamed) #colnames
    
    # find intersection
    AD_intersect = list(set(AD_GOLD) & set(AD_NAMES))
    DB_intersect = list(set(DB_GOLD) & set(DB_NAMES))
    mix_index = [0,1,2,3]
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
        dicts = {"is_{}".format(i): {} for i in mix_index}
        
        all_pairs = list(itertools.product(AD_intersect, DB_intersect))
        print(datetime.datetime.now()) 
        for pair in all_pairs:
            # 0 is AD ; 1 is DB
            # get this pair from IS_normed
            
            s = IS_normed.loc[[pair[0]+"_BC-1",pair[0]+"_BC-2"], [pair[1]+"_BC-1", pair[1]+"_BC-2"]]
            s = s.unstack().reset_index()
            s.columns = ["DB","AD", "Score"]
            scores = s.sort_values(by="Score", ascending=False).reset_index(drop=True)
            
            for i in mix_index:
                d_name = "is_{}".format(i)
                try:
                    s = scores.loc[i,:]
                except Exception:
                    continue    
            
                if str(s.Score) != "nan":
                    dicts[d_name][(s.AD, s.DB)] = s.Score
        
        print(datetime.datetime.now())
        AD_freq = AD_freq.to_frame()
        DB_freq = DB_freq.to_frame()
        
        AD_freq["AD_name"] = AD_freq.index
        DB_freq["DB_name"] = DB_freq.index
        for i in dicts.keys(): # calculate MCC for each index
        
            # convert i to df 
            # (AD, DB): score 
            sort_is = pd.Series(dicts[i]).reset_index()
            sort_is.columns = ["AD_name", "DB_name", "Score"]
            
            # convert values to numbers
            sort_is[["Score"]] = sort_is[["Score"]].apply(pd.to_numeric)
        
            # split bc and orf name
            sort_is[["AD","AD_BC"]] = sort_is["AD_name"].str.split("_",expand=True)
            sort_is[["DB","DB_BC"]] = sort_is["DB_name"].str.split("_",expand=True)
        
            # merge for interactions
            sort_is["Interactions"] = sort_is[["AD", "DB"]].apply(lambda x: "-".join(x),axis=1)

            # how many interactions we have?
            print sort_is.shape
            
            # compare gold with our set
            sort_is["ishit_interaction"] = sort_is.Interactions.isin(gold_st.Interactions).astype(int)
            
            # print hit.shape 
            # filter the sort_is
            # AD has to appear in AD_GOLD at least once
            # union of AD in our set & AD in GOLD 
            # DB has to appear in DB_GOLD at least once
            # --> both col has to be in the union set
        
            # screen = sort_is.loc[sort_is.AD.isin(AD_GOLD) | sort_is.DB.isin(DB_GOLD)]
            
            sort_is["is_AD"] = sort_is.AD.isin(AD_GOLD)
            sort_is["is_DB"] = sort_is.DB.isin(DB_GOLD)
            
            sort_is = pd.merge(sort_is, AD_freq, on="AD_name")
            sort_is = pd.merge(sort_is, DB_freq, on="DB_name")

            sort_is.rename(columns={"0_x":"AD_mfreq", "0_y":"DB_mfreq"}, inplace=True)

            #sort_is = screen.loc[(screen.AD_mfreq > row_cut) & (screen.DB_mfreq >col_cut)]
            
            sort_is["AD_cut"] = sort_is.AD_mfreq > row_cut
            sort_is["DB_cut"] = sort_is.DB_mfreq >col_cut
            
            sort_is["screen"] = (sort_is.is_AD & sort_is.is_DB & sort_is.AD_cut & sort_is.DB_cut).astype(int)

            sort_is = sort_is.sort_values(by="Score", ascending=False)  
            sort_is = sort_is.reset_index(drop=True)
            
            # hit set 
            y_pred = sort_is.ishit_interaction.tolist()
            # screen set
            y_lab = sort_is.screen.tolist()

            # scores
            y_score = sort_is.Score.tolist()
            
            
            # get confusion matrix
            cnf_matrix = confusion_matrix(sort_is.screen.tolist(), sort_is.ishit_interaction.tolist())
            print sum(sort_is.ishit_interaction.tolist())
            #print sort_is.shape
            print cnf_matrix
            tn, fp, fn, tp = cnf_matrix.ravel()
            print tn, fp, fn, tp
            
            precision, recall, _ = precision_recall_curve(y_pred, y_score)
            mcc = matthews_corrcoef(y_lab, y_pred)
            print precision
            print recall
            print mcc
            output_file = "/home/rothlab/rli/www/html/test.png"
            plot_prc(precision, recall, output_file)

            #print(datetime.datetime.now())
            #print screen
            #print screen.shape
            # from pre matrix, select those pairs 
            # build a screen set
            # If presort score sucks for this pair?
            # compare to gold standard
            # find intersections
            break
        break

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='BFG-Y2H Scores')
    parser.add_argument("--counts", help="Path for read counts files")
    parser.add_argument("--output", help="Path for output files")

    args = parser.parse_args()
    output = args.output
    counts = args.counts
    
    print "Concatenate samples ... "
    
    samples = {}

    for d in os.listdir(counts):
        sample_name = d.split("_")[0]
        matrix = d+"/combined_counts.csv"
        try:
            samples[sample_name].append(matrix)
        except Exception:
            samples[sample_name] = [matrix]
    print len(samples.keys())

    weights = np.arange(0, 2.4, 0.4)
    floor_perc = np.arange(5,25,5)
    #opt_mix = [0,1,2,3]

    os.chdir(counts)
    yi = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/YI_1.txt"
    yi_GOLD = load_YI1(yi)
    for key in samples.keys():
        for matrix in samples[key]:
            fname = "./"+matrix
            if "pre" in matrix:
                GFP_pre = pd.read_table(fname, sep=",", index_col=0)
            elif "med" in matrix:
                GFP_med = pd.read_table(fname, sep =",", index_col=0)
            elif "high" in matrix:
                GFP_high = pd.read_table(fname, sep =",", index_col=0)
        print GFP_pre.shape, GFP_med.shape, GFP_high.shape
        
        score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, yi, output)
        
        break
    #score_main(weights, floor_perc,opt_mix)
    


