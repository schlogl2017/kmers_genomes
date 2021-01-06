# usr/bin/env python
# Gets some basic stats from a fasta file (can be compressed, '.gz).
# Returns a dictionary with the name of the file and its statistics: name, number of data(seqs), length of sequence,
# base composition, N counts, CG content, alternative bases count
import sys
import fasta_parser
from sequence_utils import base_stats, get_gc_content
from system_utils import get_fasta_files

if len(sys.argv) < 3:
    print('The argument number is insufficient to correct run the script')
    print('USAGE: python fasta_stats.py <dir_name> <alphabet>')
    sys.exit()


def fasta_stats_genomes(filename):
    names = filename.split('/')[-1]
    alt_bases = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']
    data_lst = []
    for name, seq in fasta_parser.parse_fasta(filename):
        seq_len = len(seq)
        bases = base_stats(seq.upper(), alphabet, False, True)
        cg_content = get_gc_content(seq.upper())
        ns = seq.upper().count('N')
        alt = {al: seq.upper().count(al) for al in alt_bases}
        num_seqs = fasta_parser.fasta_item_counter(filename)
        headers = fasta_parser.get_fasta_headers(filename)
        data_lst.append((names, seq_len, num_seqs, bases, ns, cg_content, alt, headers))
    return data_lst


dir_name = sys.argv[1]
alphabet = sys.argv[2]
filenames = get_fasta_files(dir_name)


for filename in filenames:
    for data in fasta_stats_genomes(filename):
        name, seq_len, seq_num, bases, ns, gc, alt, h = data
        print(f'Name: {name}')
        print(f'Length: {seq_len}')
        print(f'Number of sequences: {seq_num}')
        print(f'Headers {h}')
        print(f'Base composition: {bases}')
        print(f'Num unknow bases (Ns): {ns}')
        print(f'CG content: {round(gc * 100, 4)} %')
        print(f'Alternative bases: {alt}\n')
