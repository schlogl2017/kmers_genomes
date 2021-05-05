#!/usr/bin/env python
# coding: utf-8
import os
import gzip
import argparse
import random
import textwrap
from fasta_parser import parse_fasta, str_punctuation_strip
from sequence_utils import get_name_fasta
from system_utils import get_fasta, make_me_a_folder, get_full_name, get_species_name


def parse_arguments():
    """Parse the command line arguments to the the count k-mer script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='A script to shuffle genomes.')
    parser.add_argument('-d',
                        '--dir',
                        metavar='dir_name',
                        type=str,
                        required=True,
                        dest='dir_name',
                        help='directory name. Eg. Data')                  
    parser.add_argument('-ssd',
                        '--sub_sub_dir',
                        type=str,
                        dest='sub_sub_dir',
                        help='subdirectory of the files. Eg. chromosomes or plasmids')
    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name to the output file. Eg. Results')
    parser.add_argument('-f',
                        '--file',
                        type=str,
                        dest='file_txt',
                        help='species list  as txt.')
    return parser.parse_args()


def main():
    cwd = os.getcwd()
    print(f'The working directory: {cwd}\n')
    start_time = time.process_time()
    opt = parse_arguments()
    dir_name = opt.dir_name
    sub_sub_dir = opt.sub_sub_dir
    dir_out = opt.dir_out
    file_txt = opt.file_txt
    
    spc_names = get_species_name(file_txt)
    filenames = list(get_fasta(spc_names, dir_name, sub_sub_dir))
    
    cnt_files = 0
    for filename in filenames:
        print(f'Working with the file: {filename}\n')
        for name, seq in parse_fasta(filename):
            name = get_name_fasta(name)
            spc = name.split('_')[2]
            sequence = [seq[i:1+2] for i in range(0,len(seq),2)]
            random.shuffle(sequence)
            sequence = ''.join(sequence)
            cnt_files += 1
            dir_save = get_full_name(dir_out, spc, sub_sub_dir)
            print(f'Saving the schuffled genome at: {dir_save}\n')
            if os.path.exists(dir_save):
                pass
            else:
                make_me_a_folder(dir_save)
            with gzip.open(f'{dir_save}/{name}_shuffled.fna.gz', 'wt') as fout:
                fout.write(f'>{name}_genome_shuffled\n')
                fout.write(textwrap.fill(seq, 80))
    end = time.process_time()
    total_time = end - start_time
    print(f'The script takes {total_time} to finish!')
    print(f'Where read and manipulated {cnt_files} files')
    print('Done!')
   
