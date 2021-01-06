#!usr/bin/env python
# -*- coding: utf-8 -*-


from collections import defaultdict, Counter
import itertools
import alphabet
from sequence_utils import get_reverse_complement, count_umbiguous_bases, get_chunks, get_sequence_from_scaffolds
from fasta_parser import fasta_item_counter, parse_fasta


def get_all_possible_kmers(alphabet, min_k, max_k):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet."""
    kmers = [''.join(letters) for n in range(min_k, max_k + 1)
             for letters in itertools.product(alphabet, repeat=n)]
    return kmers


def count_kmers(sequence, alphabet, min_k, max_k):
    alphabet = set(alphabet)
    counts = defaultdict(int)
    for kmer in get_kmers_from_sequence(sequence, min_k, max_k):
        if set(kmer).issubset(alphabet):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def get_kmer_counts(kmer_list, kmer_counts):
    counts = defaultdict(int)
    for kmer in kmer_list:
        counts[kmer] = counts.get(kmer, 0) + kmer_counts[kmer]
    return counts


def get_kmers_from_sequence(sequence, min_k, max_k):
    """
    Generate all DNA k-mers over the entirety of a sequence.
    Inputs:
    sequence - string where all kmers will be checked
    min_k: minimum DNA kmer length (int)
    max_k: maximum DNA kmer length (int)
    Output:
    yields all DNA kmers (str) of length min_k to max_k
    """
    limits = range(min_k, max_k + 1)
    for i in range(0, len(sequence) - max_k + 1):
        for j in limits:
            yield sequence[i:i+j]


def kmer_positions(sequence, alphabet, k):
    """ returns the position of all k-mers in sequence as a dictionary"""
    mer_position = defaultdict(list)
    for i in range(1, len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if all(base in set(alphabet) for base in kmer):
            mer_position[kmer] = mer_position.get(kmer, []) + [i]
    # combine kmers with their reverse complements
    pair_position = defaultdict(list)
    for kmer, pos in mer_position.items():
        krev = get_reverse_complement(kmer)
        if kmer < krev:
            pair_position[kmer] = sorted(pos + mer_position.get(krev, []))
        elif krev < kmer:
            pair_position[krev] = sorted(mer_position.get(krev, []) + pos)
        else:
            pair_position[kmer] = pos
    return pair_position


def get_kmer_count_slide_window(sequence, alphabet, window, step, min_k, max_k):
    slide_mer_count = defaultdict(Counter)
    for chunk, s, e in get_chunks(sequence, window, step):
        pos = '_'.join((str(s), str(e)))
        slide_mer_count[pos].update(count_kmers(chunk, alphabet, min_k, max_k))
    return slide_mer_count


def get_kmer_clumps(sequence, alphabet, k, window, times):
    clumps = defaultdict(list)
    kmers = kmer_positions(sequence, alphabet, k)
    for kmer, pos in kmers.items():
        clumps[kmer] = clumps.get(kmer, [])
        for i in range(len(pos) - times):
            end = i + times -1
            while (pos[end] - pos[i]) <= window - k:
                end += 1
                if end >= len(pos):
                    break
            if end - i >= times:
                clumps[kmer].append((pos[i], end - i))
    return clumps


if __name__ == '__main__':
    alphabet = alphabet.iupac_dna
    filename = "Data/test/Abditibacterium/GCF_002973605.1_ASM297360v1_genomic.fna.gz"
    min_k, max_k = 4, 4
    kmers = []
    if fasta_item_counter(filename) > 1:
        sequence = get_sequence_from_scaffolds(filename)
        for kmer in get_kmers_from_sequence(sequence, min_k, max_k):
            kmers.append(kmer)
    else:
        for name, seq in parse_fasta(filename):
            sequence = seq
        for kmer in get_kmers_from_sequence(sequence, min_k, max_k):
            kmers.append(kmer)
    print(f'Sequence length: {len(sequence)}')
    print(f'Unknow base counts: {count_umbiguous_bases(sequence)}')
    print(f'Theoretical number of kmers of length 4: {4 ** max_k}')
    print(f'Number of Kmer founds in the sequence: {len(set(kmers))}')
    print(set(kmers))
