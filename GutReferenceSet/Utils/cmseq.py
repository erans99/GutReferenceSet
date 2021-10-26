import os
import pandas as pd


def read_cmseq(path):
    df = pd.read_csv(path, index_col=0)
    df.index = [os.path.basename(i).replace('.fa', '') for i in df.index]
    return df
