import os
import pandas as pd

class ParseSam(object):

    def __init__(self, sam_file):
        self._sam_file = sam_file

    def _Parse(self, self._sam_file):
        """ 
        Read the header of sam file and save it in a separate file
        Read the content of sam file and load as a dictionary
        Return loaded dictionary
        """
        header = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", 
                "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"]

        content = []
        filename = os.path.basename(self._sam_file).split(".")[0]
        with open(self._sam_file, "r") as f, open(filename+".header", "r") as h:
            for line in f:
                if "@" in line:
                    h.write(line+"\n")
                else:
                    content.append(line)

        file_content = pd.DataFrame(nt, columns=headers)
        return file_content 
    def _Stat(self):
        """
        Return the stat of the sam file
        """
        # total samples identified
        pass

    def _FilterQNAME(self, self._file_content, QNAME=""):
        """
        """
        pass     
