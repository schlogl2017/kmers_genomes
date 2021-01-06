!#usr/bin/env python
import os
import glob
import pandas as pd
from functools import reduce


pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


extension = 'csv'


dfs = [
    pd.read_csv(csv_file, sep=',') for csv_file in
    [f for f in glob.glob(f'Data/test/*.{extension}')]
]
merged = reduce(lambda left, right: pd.merge(left, right, on='kmer'), dfs)
merged.to_csv('out.csv', sep=',', index=False)
