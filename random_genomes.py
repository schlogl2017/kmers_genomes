#!/usr/bin/env python
# coding: utf-8


import os
import sys
import gzip
import time
import random
import argparse
from fasta_parser import parse_fasta, str_punctuation_strip
from system_utils import get_files, make_me_a_folder


def write_shuffled_genome(filename, shuffled_filename):
    """Write a shuffled version of the sequence to a FASTA file."""
    for name, sequence in parse_fasta(filename):
        name = '_'.join(str_punctuation_strip(name)[1:4])
        seq = sequence.upper()
        seq_list = list(seq)
        random.shuffle(seq_list)
        shuffled_seq = ''.join(seq_list)
        with gzip.open(shuffled_filename, 'wt') as fout:
            fout.write(f'>shuffled {name} sequence\n')
            fout.write(shuffled_seq)


def parse_arguments():
    """Parse the command line arguments to the the count k-mer script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to shuffle genomes.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Path to the files')
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
    for filename in filenames:
        write_shuffled_genome(filename, filename[:-7] + '_ranomized.fa.gz')
        cnt_files += 1
        print(f'Working with the file: {filename}')
    end = time.process_time()
    total_time = end - start_time
    print(f'The script takes {total_time} to finish!')
    print(f'Where read and manipulated {cnt_files} files')
    print('Done!')


if __name__ == "__main__":
    sys.exit(main())




