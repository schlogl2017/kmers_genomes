#!usr/bin/env python
#  python kmer_count_genomes.py -p Data/test -d Results/bacteria -o bacterias -ka 5 -ki 5
#


import os
import sys
import time
import argparse
import pandas as pd
from collections import defaultdict
import fasta_parser
from get_kmers import get_all_possible_kmers, count_kmers, get_kmer_counts
from system_utils import get_files, make_me_a_folder


def parse_arguments():
    """Parse the command line arguments to the the count k-mer script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to count kmer in genomes.')
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
                        default='ACGT',
                        dest='alphabet',
                        help='The allowed characters tha compound the sequence.')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        dest='outfile',
                        help='Name for output files.')
    parser.add_argument('-ka',
                        '--max_k',
                        type=int,
                        default=2,
                        action='store',
                        dest='max_k',
                        help='Specify the maximum kmer/palindrome length.')
    parser.add_argument('-ki',
                        '--min_k',
                        type=int,
                        default=2,
                        action='store',
                        dest='min_k',
                        help='Specify the minimum kmer/palindrome length.')
    return parser.parse_args()


def main():
    """Parses options from the command line.
    Computes the k-mers to test (either palindromes or all k-mers).
    Computes the counts of k-mers in fasta files, and add the reverse complements
    of the sequence data to the counts.
    Computes the kmers/palindromes statistics (expected value, z-scores and e-values),
    And if definide by user prints the results to stdout, else save to a csv file.
    """
    cwd = os.getcwd()
    print(f'The working directory: {cwd}')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.path
    filenames = get_files(dir_name)
    outfile = opt.outfile
    dir_out = opt.dir_out
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)
    cnt_files = 0
    dict_counts = defaultdict(dict)
    for filename in filenames:
        print(f"Working on the file: {filename}")
        file_n = filename.split('/')[-1]
        for name, seq in fasta_parser.parse_fasta(filename):
            dict_counts[file_n].get(file_n, {})
            seq = seq.upper()
            kmer_list = get_all_possible_kmers(opt.alphabet, opt.min_k, opt.max_k)
            kmer_counts = count_kmers(seq, opt.alphabet, opt.min_k, opt.max_k)
            counts = get_kmer_counts(kmer_list, kmer_counts)
            dict_counts[file_n].update(counts)
            df = pd.DataFrame(dict_counts)
            df.to_csv(f'{dir_out}/{outfile}_kmers{opt.min_k}_to_{opt.max_k}.csv')
        cnt_files += 1
    end = time.process_time()
    total_time = end - start_time
    print(f'The script takes {total_time} to finish!')
    print(f'Where read and manipulated {cnt_files} files')
    print('Done!')


if __name__ == "__main__":
    sys.exit(main())
