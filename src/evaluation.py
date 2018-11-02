import score
import noz_score
import param
import pandas as pd
import numpy as np

def load_litbm(lit):
    """
    load litbm to pandas dataframe
    input cols: AD DB
    output cols: AD DB Interactions
    """
    litbm = pd.read_table(lit)
    litbm["Interactions"] = litbm[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return litbm

def noz_main(litbm, score_matrix, max_rank, max_rho):
    # get normalized score by using max_pho
    norm_s = noz_score.norm_score(score_matrix, max_pho)
    # get optimized rank by using max_rank
    r = int(max_rank.split("_")[1])
    output_rank = noz_score.getrank(norms, r)
    mcc = noz_score.get_screen(output_rank, litbm)
    print mcc

def dk_main(litbm, max_weight, max_rank, max_floor):
    pass


if __name__ == "__main__":
    
    test_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/rerun_analysis/yAD1DB1/"
    litbm=load_litbm(param.litBM13)
    print litbm
    #noz_main()

