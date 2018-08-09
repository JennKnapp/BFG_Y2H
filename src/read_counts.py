import numpy as np
import pandas as pd
import os
import datetime
from param import *
#from main import *
import logging
import logging.config

#logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
analysis_log = logging.getLogger("analysis")


class Read_Count(object):

    def __init__(self, AD_GENES, DB_GENES, r1, r2):
        self._r1 = r1 # sorted sam file for r1
        self._r2 = r2 # sorted sam file for r2
        
        self._ad_genes = AD_GENES # list of AD gene names
        self._db_genes = DB_GENES # list of DB gene names 
        
        ####### log ########
        analysis_log.info("r1 sam file: %s", r1)
        analysis_log.info("r2 sam file: %s", r2)
        
    def _BuildMatrix(self):
        """
        Build empty up tag matrix with col = DB_size ; row = AD_size 
        Build empty dn tag matrix with col = DB_size ; row = AD_size
        """
        self._uptag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
        self._dntag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
        print "Matrix build" 
        analysis_log.info("Matrix build with AD: %s, DB: %s", str(self._uptag_matrix.shape[0]), str(self._uptag_matrix.shape[1]))

    def _ReadCounts(self):
        
        f1 = open(self._r1, "rb")
        f2 = open(self._r2, "rb")
        
        i=True
        lines = 0
        print datetime.datetime.now()
        while 1:
            lines+=1
            if lines % 1000000 ==0:
                print lines
                print datetime.datetime.now()
            
            r1_line = f1.readline() # AD
            r2_line = f2.readline() # DB
            
            if r1_line == ""  or r2_line == "":
                i = False
                print "end of file"
                break

            r1_line = r1_line.strip().split("\t")
            r2_line = r2_line.strip().split("\t")
            
            # both files are sorted by name, if name is different, log error
            if r1_line[0] != r2_line[0]:
                #log error and exit
                i = False
                analysis_log.error("# READ ID DOES NOT MATCH #")
                analysis_log.error("From read one (AD): %s", r1_line[0])
                analysis_log.error("From read two (DB): %s", r2_line[0])
                break

            if int(r1_line[4]) < 3 or int(r2_line[4]) < 3: # check quality
                continue

            if r1_line[2] == "*" or r2_line[2] =="*": # if one of the read didnt map
                continue
            
            r1_name = r1_line[2].split(";")
            r2_name = r2_line[2].split(";")

            if r1_name[-1] == r2_name[-1]:
                #print r1_name
                if r1_name[-1] == "dn":
                    self._dntag_matrix.loc[r1_name[1], r2_name[1]] +=1
                else:
                    self._uptag_matrix.loc[r1_name[1], r2_name[1]] +=1
        
        print "file read"
        f1.close()
        f2.close()
        uptag_file = "./uptag_rawcounts.csv"
        dntag_file = "./dntag_rawcounts.csv"
        self._dntag_matrix.to_csv(dntag_file)
        self._uptag_matrix.to_csv(uptag_file)
        print "matrix saved"
        #return self._dntag_matrix, self._uptag_matrix


def RCmain(r1, r2, AD_genes, DB_genes):
    analysis_log.info("test")
    rc = Read_Count(AD_genes, DB_genes, r1, r2)
    # create empty matrix
    rc._BuildMatrix()
    
    # create 
    rc._ReadCounts()

if __name__ == "__main__":
    # test rc main
    r1_csv = "/home/rothlab/rli/02_dev/08_bfg_y2h/yAD1DB4_output/yAD1DB4_GFP_high/yAD1DB4_GFP_high_R1_AD_BC_noh.csv"
    r2_csv = "/home/rothlab/rli/02_dev/08_bfg_y2h/yAD1DB4_output/yAD1DB4_GFP_high/yAD1DB4_GFP_high_R2_DB_BC_noh.csv"
    AD_genes, DB_genes = read_summary(AD_summary, DB_summary, AD_group="G1", DB_group="G4")
    print AD_genes
    print DB_genes
    dn_matrix, up_matrix = RCmain(r1_csv, r2_csv, AD_genes, DB_genes)
