#!usr/bin/env python
# -*- coding: utf-8 -*-
#!usr/bin/env python
# python run_parallel.py -p Data/bacteria -d Results -k 2

import time
import gzip
import os
import sys
import argparse
import pandas as pd
from joblib import Parallel, delayed
from system_utils import make_me_a_folder, get_dir_name, get_fasta_files


def indexfasta(filename):
    # Open file
    infile = gzip.open(filename)

    # Set chunksize
    chunksize = 1024 * 1024
    filepos = 0
    headstart = list()
    headend = list()

    # Read chunck size
    while True:
        content = infile.read(chunksize)
        # print(content)
        if len(content) == 0:
            break

        # Find Header Start
        chunkpos = 0
        while chunkpos != -1:
            chunkpos = content.find(b'>', chunkpos)
            if chunkpos != -1:
                headstart.append(chunkpos + filepos)
                chunkpos += 1

        # Find Header End
        for i in range(len(headend), len(headstart)):
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()

    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart) - 1, 0, -1):
        if headend[i] == headend[i - 1]:
            del headstart[i]
            del headend[i]
    headstart.append(filepos)
    fastaindex = list()
    for i in range(len(headend)):
        seq_start = headend[i] + 1
        seq_end = headstart[i + 1] - 1
        fastaindex.append((seq_start, seq_end, seq_end - seq_start))
    return fastaindex


def indexsequence(seq):
    seq = seq.upper()
    pointer = 0
    seqindex = list()
    while len(seq) > pointer:
        # Find start of seq
        potenstart = [seq.find(b'A', pointer),
                      seq.find(b'T', pointer),
                      seq.find(b'C', pointer),
                      seq.find(b'G', pointer)]
        realstart = min(potenstart)
        if realstart == -1:
            # happens rarely, so slow code is ok
            potenstart = [i for i in potenstart if i > -1]
            if len(potenstart) == 0:
                break
            realstart = min(potenstart)
        realend = seq.find(b'N', realstart)
        if realend == -1:
            realend = len(seq)
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex


def find_kmers(fasta, idx, transtable, kmer_len):
    # Read sequence
    infile = gzip.open(fasta, 'rb')
    infile.seek(idx[0])
    seq = infile.read(idx[1] - idx[0] + 1).translate(transtable, b'\r\n\t ')
    infile.close()
    subdict = dict()
    # Index sequence
    seqindex = indexsequence(seq)
    # Slide through sequences and add to dictionary
    for start, stop in seqindex:
        for i in range(start, stop - kmer_len + 1):
            kmer = seq[i:i + kmer_len]
            if kmer not in subdict:
                subdict[kmer] = 1
            else:
                subdict[kmer] += 1
    return subdict


def parse_arguments():
    """Parse the command line arguments to the the count k-mer script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. The resulting results of count kmers is a csv file with all
    kmers and it's counts.
    """
    parser = argparse.ArgumentParser(description='A script to count kmer in genomes.')
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
    parser.add_argument('-k',
                        '--length',
                        type=int,
                        default=2,
                        action='store',
                        dest='k',
                        help='Specify the kmer length.')
    return parser.parse_args()


def main():
    cwd = os.getcwd()
    print(f'The working directory: {cwd}')
    opt = parse_arguments()
    dir_name = opt.path
    filenames = get_fasta_files(dir_name)
    #outfile = opt.outfile
    dir_out = opt.dir_out
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)
    kmer_len = opt.k
    transtable = bytes.maketrans(b'ACGTMRYKVHDBWmrykvhdbxnsw', b'ACGTNNNNNNNNNNNNNNNNNNNNN')
    # n_worker = os.cpu_count()
    # INDEXING
    start = time.time()
    cnt_files = 0
    for file in filenames:
        path_out = os.path.join(*get_dir_name(file))
        name = file.split('/')[2]
        indexes = indexfasta(file)
        index_time = time.time() - start
        # FIND KMER
        start = time.time()
        results = Parallel(n_jobs=4)(delayed(find_kmers)(file, z, transtable, kmer_len) for z in indexes)
        search_time = time.time() - start
        # MERGE DICTS
        final_dict = dict()
        start = time.time()
        for r in results:
            for mer in r.keys():
                if mer in final_dict:
                    final_dict[mer] += r[mer]
                else:
                    final_dict[mer] = r[mer]
            cnt_files += 1
        df = pd.DataFrame([(x.decode("utf-8"), cnt) for x, cnt in final_dict.items()])
        fullname = os.path.join(dir_out, name, f'kmer{kmer_len}')
        if not os.path.exists(fullname):
            os.makedirs(fullname)
        df.to_csv(f'{fullname}/{name}_kmer{kmer_len}.csv',
                  sep=',',
                  header=['kmer', 'counts'],
                  index=False)
        dict_sum_time = time.time() - start
        del results
        del r
        print(f"<b> Results for file <b>: {name}")
        print("Total time: ", round(index_time + search_time + dict_sum_time, 3))
        print(f"Total number of files: {cnt_files}")
    print('Done')


if __name__ == "__main__":
    sys.exit(main())


