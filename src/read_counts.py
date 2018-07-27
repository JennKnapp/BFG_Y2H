import numpy as np
import pandas as pd
from parse_sam import *
from param import *

class Read_Count(object):

    def __init__(self, AD_GENES, DB_GENES, r1, r2):
        self._r1 = r1 # sorted sam file for r1
        self._r2 = r2 # sorted sam file for r2
        self._ad_size = AD_SIZE # number of genes in AD group
        self._db_size = DB_SIZE # number of genes in DB group
        
        self._ad_genes = AD_GENES # list of AD gene names
        self._db_genes = DB_GENES # list of DB gene names 

    def _BuildMatrix(self):
        """
        Build empty up tag matrix with col = DB_size ; row = AD_size 
        Build empty dn tag matrix with col = DB_size ; row = AD_size
        """
        self._uptag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
        self._dntag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
                

    def _ReadCounts(self):
        
        f1 = open(self._r1, "rb")
        f2 = open(self._r2, "rb")
        
        while 1:
            r1_line = f1.readline() # AD
            r2_line = f2.readline() # DB
            
            if "@HD" or "@PG" in r1_line:
                continue
            
            if not r1_line or not r2_line:
                break

            r1_line = r1_line.split("\t")
            r2_line = r2_line.split("\t")
            
            # both files are sorted by name, if name is different, log error
            if r1_line[0] != r2_line[0]:
                #log error and exit
                break

            if int(r1_line[2]) < 3 or int(r2_line[2]) < 3: # check quality
                continue

            if r1_line[4] == "*" or r2_line[4] =="*": # if one of the read didnt map
                continue
            
            r1_name = r1_line[4].split(";")
            r2_name = r2_line[4].split(";")

            if r1_name[-1] == r2_name[-1]:
                if r1_name[-1] == "dn":
                    self._dntag_matrix.loc[[r1_name[1]], [r2_name[1]]] +=1
                else:
                    self._uptag_matrix.loc[[r1_name[1]], [r2_name[1]]] +=1


def RCmain(r1, r2, AD_genes, DB_genes):

    rc = Read_Count(AD_genes, DB_genes, r1, r2)
    
    # create empty matrix
    empty_matrix = rc._BuildMatrix()
    # create 
    count_reads = rc._ReadCounts()
