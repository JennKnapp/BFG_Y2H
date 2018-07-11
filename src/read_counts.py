import numpy as np

def parse_sam_file(sam_file, AD_size, DB_size):

    """
    Count barcodes in sam_file
    """
    uptag_matrix = np.zeros(AD_size, DB_size)
    dntag_matrix = np.zeros(AD_size, DB_size)

    with open(sam_file, "r") as f:
        for line in f:
            if "@SQ" or "@HD" in line: continue


    pass
