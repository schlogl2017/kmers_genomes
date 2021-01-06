#!/usr/bin/env python
# coding: utf-8
# python kmer_palindromes_analysis.py -p Data/test -d Results -pa -ka 6 -ki 6
# python kmer_palindromes_analysis.py -p Data/test -d Results -km -ka 6 -ki 6
# python kmer_palindromes_analysis.py -p Data/test -d Results  -a -ka 6 -ki 6
#  python kmer_palindromes_analysis.py -p Data/test -d Results -sl -w 20000 -s 10000 -ka 6 -ki 6

import os
import sys
import time
import argparse
import fasta_parser
from system_utils import get_files, make_me_a_folder, get_dir_name, get_full_name
from get_kmers import get_all_possible_kmers
from kmer_palindromes_count import count_kmers, kmers_frequencies, get_kmer_count_slide_window, print_results_stats
from kmer_palindrome_stats import *
from sequence_utils import count_umbiguous_bases, get_reverse_complement
from kmer_palindrome_save_data import *
from get_palindromes import get_palindromes
from alphabet import iupac_dna


def parse_arguments():
    """Parse the command line arguments to the genome_palindromes_analysis script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to analyze palindromes and kmer in genomes.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Path to the files')
    parser.add_argument('-alp',
                        '--alph',
                        type=str,
                        default=iupac_dna,
                        dest='alphabet',
                        help='The allowed characters tha compound the sequence.')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        dest='output',
                        help='name output files.')
    parser.add_argument('-e',
                        '--max_e',
                        type=float,
                        default=0.01,
                        action='store',
                        dest='max_e',
                        help='Threshold e value for report over/under represented palindromes.')
    parser.add_argument('-ka',
                        '--max_k',
                        type=int,
                        default=6,
                        action='store',
                        dest='max_k',
                        help='Specify the maximum kmer/palindrome length.')
    parser.add_argument('-ki',
                        '--min_k',
                        type=int,
                        default=6,
                        action='store',
                        dest='min_k',
                        help='Specify the minimum kmer/palindrome length.')
    parser.add_argument('-w',
                        '--window',
                        type=int,
                        action="store",
                        dest='window',
                        help='Window size for the sequences chunks.')
    parser.add_argument('-s',
                        '--step',
                        type=int,
                        action="store",
                        dest='step',
                        help='Step size for the sequences chunks.')
    parser.add_argument('-sl',
                        '--slide',
                        action="store_true",
                        dest='slide',
                        help='Count a pattern in a sequence.')
    parser.add_argument('-a',
                        '--all',
                        action='store_true',
                        dest='all',
                        help='Test all possible k-mers (min_k <= k <= max_k.')
    parser.add_argument('-km',
                        '--kmer',
                        action='store_true',
                        dest='kmer',
                        help='Test all possible k-mers (min_k <= k <= max_k.')
    parser.add_argument('-pa',
                        '--pal',
                        action='store_true',
                        dest='pal',
                        help='Test all possible palindromes (min_k <= k <= max_k.')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    Computes the k-mers to test (either palindromes or all k-mers).
    Computes the counts of k-mers in fasta files, and add the reverse complements
    of the sequence data to the counts.
    Computes the k-mers/palindromes statistics (expected value, z-scores and e-values),
    And if definide by user prints the results to stdout, else save to a csv file.
    """
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.path
    filenames = get_files(dir_name)
    outfile = opt.output
    dir_out = opt.dir_out
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)
    cnt_files = 0
    for filename in filenames:
        print(f"Working on the file: {filename}")
        file_n = filename.split('/')[-1]
        sub, subsub = get_dir_name(filename)
        for name, seq in fasta_parser.parse_fasta(filename):
            seq = seq
            len_seq = len(seq) - count_umbiguous_bases(seq)
            if opt.kmer:
                kmer_counts = count_kmers(seq,
                                          opt.alphabet,
                                          opt.min_k - 2,
                                          opt.max_k)
                kmer_list = get_all_possible_kmers(opt.alphabet,
                                                   opt.min_k,
                                                   opt.max_k)
                kmer_freqs = kmers_frequencies(kmer_counts)
                kmer_expected = get_expected_values(kmer_list,
                                                    kmer_counts)
                kmer_zscores = get_z_scores(kmer_list,
                                            kmer_counts,
                                            kmer_expected,
                                            len_seq)
                kmer_pvalues = get_pvalues(kmer_list,
                                           kmer_zscores)
                kmer_evalues = get_evalues(kmer_list,
                                           kmer_pvalues)
                kmer_scores = get_scores(kmer_list,
                                         kmer_counts,
                                         kmer_expected)
                kmer_nscores = get_new_scores(kmer_list,
                                              kmer_counts,
                                              kmer_expected)
                kmer_odds_ratio = get_odds_ratio(kmer_list,
                                                 kmer_freqs)
                kmer_diff = get_difference(kmer_list,
                                           kmer_counts,
                                           kmer_expected)
                kmer_lod = get_log_odds(kmer_list,
                                        kmer_counts,
                                        kmer_expected)
                kmer_data = get_kmer_statistics(kmer_list,
                                                kmer_counts,
                                                kmer_expected,
                                                kmer_zscores,
                                                kmer_evalues,
                                                kmer_odds_ratio,
                                                kmer_diff,
                                                kmer_scores,
                                                kmer_nscores,
                                                kmer_lod)
                print_results_stats(file_n,
                                    kmer_list,
                                    len_seq,
                                    opt.min_k,
                                    opt.max_k,
                                    opt.max_e,
                                    kmer_data)
                df = pd.DataFrame(kmer_data, columns=["kmer",
                                                      "Observed",
                                                      "Expected",
                                                      "Z_score",
                                                      "Evalues",
                                                      "Odds",
                                                      "Diff",
                                                      "Scores",
                                                      "NScores",
                                                      "Log_odds"])
                fullname = get_full_name(dir_out, sub, subsub)
                if not os.path.exists(fullname):
                    os.makedirs(fullname)
                df.to_csv(f"{fullname}/{file_n}_{opt.max_k}_all_kmer_stats.csv")
                with open(f"{fullname}/{file_n}_{opt.max_k}_kmer_counts.csv", 'w') as fout:
                    fout.write('Kmer,Counts\n')
                    for kmer, count in kmer_counts.items():
                        fout.write(kmer + "," + str(count) + "\n")

            if opt.pal:
                n = len_seq
                pal_list = list(get_palindromes(opt.alphabet,
                                                opt.min_k,
                                                opt.max_k))
                # counts = counts of the kmers/palindromes with min_k-2 <= k <= max_k
                pal_counts = count_kmers(seq,
                                         opt.alphabet,
                                         opt.min_k - 2,
                                         opt.max_k)
                # as palindromes are the need to count both strands
                rev_strand_cnt = dict((get_reverse_complement(kmer), cnt) for kmer, cnt in pal_counts.items())
                for kmer, cnt in rev_strand_cnt.items():
                    pal_counts[kmer] += cnt
                n *= 2
                pal_freqs = kmers_frequencies(pal_counts)
                pal_expected = get_expected_values(pal_list,
                                                   pal_counts)
                pal_zscores = get_z_scores(pal_list, pal_counts,
                                           pal_expected,
                                           len_seq)
                pal_pvalues = get_pvalues(pal_list,
                                          pal_zscores)
                pal_evalues = get_evalues(pal_list,
                                          pal_pvalues)
                pal_scores = get_scores(pal_list,
                                        pal_counts,
                                        pal_expected)
                pal_nscores = get_new_scores(pal_list,
                                             pal_counts,
                                             pal_expected)
                pal_odds_ratio = get_odds_ratio(pal_list,
                                                pal_freqs)
                pal_diff = get_difference(pal_list,
                                          pal_counts,
                                          pal_expected)
                pal_lod = get_log_odds(pal_list,
                                       pal_counts,
                                       pal_expected)
                pal_data = get_kmer_statistics(pal_list,
                                               pal_counts,
                                               pal_expected,
                                               pal_zscores,
                                               pal_evalues,
                                               pal_odds_ratio,
                                               pal_diff,
                                               pal_scores,
                                               pal_nscores,
                                               pal_lod)
                print_results_stats(file_n,
                                    pal_list,
                                    len_seq,
                                    opt.min_k,
                                    opt.max_k,
                                    opt.max_e,
                                    pal_data)
                df = pd.DataFrame(pal_data, columns=["pal",
                                                     "Observed",
                                                     "Expected",
                                                     "Z_score",
                                                     "Evalues",
                                                     "Odds",
                                                     "Diff",
                                                     "Scores",
                                                     "NScores",
                                                     "Log_odds"])

                fullname = get_full_name(dir_out, sub, subsub)
                if not os.path.exists(fullname):
                    os.makedirs(fullname)
                df.to_csv(f"{fullname}/{file_n}_{opt.max_k}_all_pal_stats.csv")
                with open(f"{fullname}/{file_n}_{opt.max_k}_palindrome_counts.csv", 'w') as fout:
                    fout.write('Palindrome,Counts\n')
                    for pal, count in pal_counts.items():
                        fout.write(pal + "," + str(count) + "\n")
            if opt.all:
                kmer_counts = count_kmers(seq,
                                          opt.alphabet,
                                          opt.min_k - 2,
                                          opt.max_k)
                kmer_list = get_all_possible_kmers(opt.alphabet,
                                                   opt.min_k,
                                                   opt.max_k)
                kmer_freqs = kmers_frequencies(kmer_counts)
                kmer_expected = get_expected_values(kmer_list,
                                                    kmer_counts)
                kmer_zscores = get_z_scores(kmer_list,
                                            kmer_counts,
                                            kmer_expected,
                                            len_seq)
                kmer_pvalues = get_pvalues(kmer_list,
                                           kmer_zscores)
                kmer_evalues = get_evalues(kmer_list,
                                           kmer_pvalues)
                kmer_scores = get_scores(kmer_list,
                                         kmer_counts,
                                         kmer_expected)
                kmer_nscores = get_new_scores(kmer_list,
                                              kmer_counts,
                                              kmer_expected)
                kmer_odds_ratio = get_odds_ratio(kmer_list,
                                                 kmer_freqs)
                kmer_diff = get_difference(kmer_list,
                                           kmer_counts,
                                           kmer_expected)
                kmer_lod = get_log_odds(kmer_list,
                                        kmer_counts,
                                        kmer_expected)
                kmer_data = get_kmer_statistics(kmer_list,
                                                kmer_counts,
                                                kmer_expected,
                                                kmer_zscores,
                                                kmer_evalues,
                                                kmer_odds_ratio,
                                                kmer_diff,
                                                kmer_scores,
                                                kmer_nscores,
                                                kmer_lod)
                get_dataframe_from_kmer_data(dir_out, outfile, opt.max_k, kmer_data)
                data_dict = defaultdict(list)
                for data in kmer_data:
                    kmer = data[0]
                    obs = data[1]
                    exp = data[2]
                    zscr = data[3]
                    evals = data[4]
                    data_dict[kmer] = data_dict.get(kmer, []) + [obs, exp, zscr, evals]
                fullname = get_full_name(dir_out, sub, subsub)
                if not os.path.exists(fullname):
                    os.makedirs(fullname)
                with open(f'{fullname}/{file_n}_all_kmers_z_scores.csv', 'w') as fout:
                    fout.write('kmer, data\n')
                    for kmer, data in data_dict.items():
                        fout.write(kmer + ',' + str(data) + '\n')

            if opt.slide:
                kmer_slide = get_kmer_count_slide_window(seq,
                                                         opt.alphabet,
                                                         opt.window,
                                                         opt.step,
                                                         opt.min_k,
                                                         opt.max_k)
                df = pd.DataFrame.from_dict(kmer_slide).fillna(0.0)
                fullname = get_full_name(dir_out, sub, subsub)
                if not os.path.exists(fullname):
                    os.makedirs(fullname)
                df.to_csv(f"{fullname}/{file_n}_slide_window.csv")
        cnt_files += 1
    end = time.process_time()
    total_time = end - start_time
    print(f'The script takes {total_time} to finish!')
    print(f'Where read and manipulated {cnt_files} files')
    print('Done!')


if __name__ == "__main__":
    sys.exit(main())
