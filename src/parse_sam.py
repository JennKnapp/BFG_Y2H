import os
import pandas as pd
from itertools import izip


class ParseSam(object):

    def __init__(self, r1_sam, r2_sam):
        self._r1_sam = r1_sam
        self._r2_sam = r2_sam

    def _LoadFile(self):
        
        pass

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
