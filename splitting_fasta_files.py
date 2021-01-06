#!usr/bin/env python
# python splitting_fasta_files.py -p Data/bacteria/Acidisarcina -d Results/test/ -o plasmids/chromosome
# -*- coding: utf-8 -*-
import time
from termcolor import colored
import os
import sys
import argparse
from fasta_parser import get_plasmids_chromosome_sequences_from_complet_fasta
from fasta_parser import write_fasta_file, str_punctuation_strip
from system_utils import make_me_a_folder, get_fasta_files


def parse_arguments():
    """Parse the command line arguments to the the split fasta script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. The resulting results of split fasta script are two
    python dictionaries (plasmids, chromosome) mapping genome/plasmids id with all
    sequences.
    """
    parser = argparse.ArgumentParser(description='A script to split compressed bacterial fasta files.')
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
    start = time.time()
    cwd = os.getcwd()
    print(colored(f'The working directory: {cwd}', attrs=['bold']))
    opt = parse_arguments()
    dir_name = opt.path
    filenames = get_fasta_files(dir_name)
    outfile = opt.outfile.split('/')
    dir_out = opt.dir_out
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)
    cnt_files = 0
    for file in filenames:
        # name of the taxon directory, ie. Acidisarcina
        name = file.split('/')[2]
        plasmids, chromosome = get_plasmids_chromosome_sequences_from_complet_fasta(file)
        plasm_names = plasmids.keys()
        full_path_plasmids = os.path.join(dir_out, name)
        if not os.path.exists(full_path_plasmids):
            os.makedirs(full_path_plasmids)
        for names in plasm_names:
            np = '_'.join(str_punctuation_strip(names))
            #print(np)
            write_fasta_file(plasmids, f'{full_path_plasmids}' + '/' + np + '.fna')
        full_path_chromosome = os.path.join(dir_out, name)
        if not os.path.exists(full_path_chromosome):
            os.makedirs(full_path_chromosome)
        chromosome_names = chromosome.keys()
        for names in chromosome_names:
            nchr = '_'.join(str_punctuation_strip(names))
            #print(nchr)
            write_fasta_file(chromosome, f'{full_path_chromosome}' + '/' + nchr + '.fna')
            cnt_files += 1
        print(colored(f"Results for file: {name}", attrs=['bold']))
        print(colored(f"Total number of files: {cnt_files}", attrs=['bold']))
    end = time.time()
    print(colored(f'Total time for the script: {end-start}', attrs=['bold']))
    print(colored('Done', attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
