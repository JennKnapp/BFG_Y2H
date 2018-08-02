import os
import pandas
import logging
import logging.config


class CalculateScore(object):

    def __init__(self, uptag_matrix, dntag_matrix):
        self._upma = uptag_matrix
        self._dnma = dntag_matrix


    def _Freq(self):
        """
        calculate frequencies for ORFs
        """
        pass


