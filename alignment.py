import pandas as pd

def build_ref(summary):
    """
    Generate fasta file from summary
    See example_summary as templete
    """

    summary = pd.read_table(summary)


