import os
import pandas as pd
from itertools import izip

class ParseSam(object):

    def __init__(self, r1_sam, r2_sam):
        self._r1_sam = r1_sam
        self._r2_sam = r2_sam

    def _LoadFile(self):
        """ 
        Read the header of sam file and save it in a separate file
        Read the content of sam file and load as a dictionary
        Return loaded dictionary
        """
        header = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", 
                "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]

        content = []
        filename = os.path.basename(self._r1_sam).split(".")[0]
        with open(self._r1_sam, "r") as f, open(self._r2_sam, "r") as f2:
            for r1, r2 in izip(f, f2):
                if "@" in r1:
                    continue
                if "@" in r2:
                    continue
                print r1.strip()
                print r2.strip()
                break
                content.append(line.strip().split("\t"))

        self._r1_content = pd.DataFrame(content, columns=header)
    
    def _ReadCounts(self):
        """
        """
        for index, row in self._r1_content.iterrows():
            pass 


    def _Stat(self):
        """
        Return the stat of the sam file
        """
        # total samples identified
        pass

    def _FilterQNAME(self, QNAME=""):
        """
        """
        pass     
