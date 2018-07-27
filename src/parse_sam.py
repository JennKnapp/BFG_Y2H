import os
import pandas as pd
from itertools import izip


class ParseSam(object):

    def __init__(self, sam_file):
        self._sam_file = sam_file

    def _LoadFile(self):
        f_content = {}
        with open(self._sam_file, "rb") as sam:
            for line in sam:
                # skip header
                if "@HD" in line: 
                    print line
                    continue
                if "@PG" in line:
                    print line
                    continue
                line = line.split("\t")
                # check qual
                if int(line[4]) < 3:
                    continue
                if line[2] == "*":
                    continue
                # save read as ID:map
                ID = ":".join(line[0].split(":")[-2:])
                name = line[2]
                f_content[ID] = name
        return f_content

    def _ReadCounts(self):
        """
        """
        pass 


    def _Stat(self):
        """
        Return the stat of the sam file
        Using samtools
        """
        # total samples identified
        pass

    def _FilterQNAME(self, QNAME=""):
        """
        """
        pass   

if __name__ == "__main__":
    test_f = "/home/rothlab/rli/02_dev/08_bfg_y2h/yAD1DB4_output/yAD1DB4_GFP_high/yAD1DB4_GFP_high_R1_AD_BC.sam"
    test = ParseSam(test_f)
    f = test._LoadFile()
    print f["15396:1039"]
