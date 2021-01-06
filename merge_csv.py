#!/usr/bin/env python
# coding: utf-8


import sys
import time
import glob
from collections import defaultdict
from functools import reduce
import numpy as np
import pandas as pd


print('Starting to process the script merge_csvs.py')
start_time = time.process_time()


def get_species_name(filenames_txt):
    species = []
    with open(filenames_txt, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


if len(sys.argv) < 3:
    print("Usage: python merge_kmer_count_csvs.py <path_to_the_files> <kmer_length> <txt_file_with_species_taxa>")
    sys.exit(1)


dir_name = sys.argv[1] if len(sys.argv) > 1 else '.'
k = sys.argv[2]
subdir = f'kmer{k}'
filename = sys.argv[3]
species = get_species_name(filename)
extension = 'csv'


spc_data = defaultdict(list)
for name in species:
    spc_data[name] = spc_data.get(name, []) + [i for
                                               i in glob.glob(f'{dir_name}/{name}/{subdir}/*.{extension}')]

names_for_df = spc_data.keys()

dfs = []
cnt = 0
for name in names_for_df:
    for filename in spc_data[name]:
        print(f'Working with the file: {filename}')
        dfs.append(pd.read_csv(filename, dtype={'kmer': str, 'count': np.int32}))
        df_final = reduce(lambda left, right: pd.merge(left, right, on='kmer'),
                          dfs).set_index('kmer').sort_index()
        df_final.to_csv(f'{dir_name}/{name}/{subdir}/{name}_kmer{k}_all_concatenate.csv', index=True)
        df_sum_rows = df_final.mean(axis=1)
        df_sum_rows.to_csv(f'{dir_name}/{name}/{subdir}/{name}_kmer{k}_concatenate.csv', index=True)
    cnt += 1

# df_merged = defaultdict(list)
# for name, df_list in dfs.items():
#     df_merged[name] = df_merged.get(name, [])
#     if len(df_list) > 1:
#         df_merged[name].append(reduce(lambda left, right: pd.merge(left,
#                                                                    right,
#                                                                    on='kmer'),
#                                       df_list).set_index('kmer').sort_index().mean(axis=1))
#     else:
#         df_merged[name].append(df_list)

# for name, data in df_merged.items():
#     name = name
#     df = pd.DataFrame(data).T
#     fullname = os.path.join(dir_name, name, subdir)
#     if not os.path.exists(fullname):
#         os.makedirs(fullname)
#     df.to_csv(f'{fullname}/{name}_kmer{k}_concatenate.csv', index=True)

end_time = time.process_time()
total_time = end_time - start_time
print(f'Total files: {cnt}')
print(f'Total time for running the script: {total_time}')
print('Script finished!')
