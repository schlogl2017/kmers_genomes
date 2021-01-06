#!usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd


def get_dict_to_csv(dir_out, outfile, k, dict_data):
    with open(f'{dir_out}/{outfile}_{k}k_counts.csv', 'w') as fout:
        fout.write('key,val' + '\n')
        for key, val in dict_data.items():
            fout.write(key + ',' + str(val) + '\n')


def get_csv_kmer_slidewindow(dir_out, fasta_name, kmer_slidewindow, k_max):
    df = pd.DataFrame(kmer_slidewindow).fillna(0.0)
    df.to_csv(f"{dir_out}/{fasta_name}_{k_max}_slide_window.csv")


def get_dataframe_from_kmer_data(dir_out, filename, k, kmerdata):
    df = pd.DataFrame(kmerdata, columns=["kmer",
                                         "Observed",
                                         "Expected",
                                         "Z_score",
                                         "Evalues",
                                         "Odds",
                                         "Diff",
                                         "Scores",
                                         "NScores",
                                         "Log_odds"])
    df.to_csv(f"{dir_out}/{filename}_{k}_all_stats.csv")
