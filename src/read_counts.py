import numpy as np
import pandas as pd
from parse_sam import *
from param import *

class Read_Count(object):

    def __init__(self, AD_GENES, DB_GENES, r1, r2):
        self._r1 = r1 # sam file for r1
        self._r2 = r2 # sam file for r2
        self._ad_size = AD_SIZE # number of genes in AD group
        self._db_size = DB_SIZE # number of genes in DB group
        self._ad_genes = AD_GENES # list of AD gene names
        self._db_genes = DB_GENES # list of DB gene names 

    def _BuildMatrix(self):
        """
        Build empty up tag matrix with col = DB_size ; row = AD_size 
        Build empty dn tag matrix with col = DB_size ; row = AD_size
        """
        uptag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
        dntag_matrix = pd.DataFrame(0, index=self._ad_genes, columns=self._db_genes)
                
        r1_sam = ParseSam(r1)
        r2_sam = ParseSam(r2)

        r1_sam_content = r1_sam._Parse()
        r2_sam_content = r2_sam._Parse()

        for index, row in r1_sam_content.iterrows():
            read_name = row["QNAME"]
            quality = row["MAPQ"]
            mapped_gene = row["RNAME"]

            # find query from r2
            # check quality (>3)
            # check if both mapped to up or both mapped to dn
            pass
