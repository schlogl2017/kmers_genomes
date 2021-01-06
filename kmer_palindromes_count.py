#!/usr/bin/env python
# coding: utf-8


import copy
from collections import Counter
import operator
import fasta_parser
from sequence_utils import base_stats, get_chunks, count_umbiguous_bases, get_reverse_complement
import alphabet
from get_palindromes import get_palindromes
from get_kmers import get_kmers_from_sequence
from kmer_palindrome_stats import *


def count_kmers(sequence, alphabet, min_k, max_k):
    counts = defaultdict(int)
    for kmer in get_kmers_from_sequence(sequence, min_k, max_k):
        if set(kmer).issubset(alphabet):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def kmers_frequencies(kmer_counts):
    freqs = defaultdict(float)
    total = sum(kmer_counts.values())
    for kmer, cnt in kmer_counts.items():
        freqs[kmer] = freqs.get(kmer, 0) + cnt / total
    return freqs


def get_kmer_probability(kmer_lst, base_freq_dict):
    """Calculates the probability of the kmer.
    Input:
    kmer - a substring of length k
    base_freq_dict - a dictionary like object with base frequencies
    Output:
    prob - probability of the kmer as a product of the bases tha
    compound the kmer string
    """
    prob = 1
    for kmer in kmer_lst:
        for base in kmer:
            prob = prob * base_freq_dict[base]
    return prob


def kmer_probabilities(kmer_list, base_freq_dict):
    probs = defaultdict(float)
    for kmer in kmer_list:
        p = get_kmer_probability(kmer, base_freq_dict)
        probs[kmer] = probs.get(kmer, 0.0) + p
    return probs


def print_results_stats(filename, mer_list, len_sequence, min_k, max_k, max_e, data_kmer):
    """
    Prints the count, expected_value z_score and e_value to stdout.
    Inputs:
    palindrome_lst - list of palindromes
    counts - a dictionary of kmer counts (key = kmer (str) and value = counts (int))
    expected - a dictionary of kmer expected values
    min_k: minimum kmer length (int)
    max_k: maximum kmer length (int)
    length_sequence = length of sequence (int)
    z_score dictionary - key as palindrome (str) and values as z_scores (float)
    e_value dictionary - key as palindrome (str) and value as e_value (float)
    Output:
    output in table format kmer observed expected Z_score E_value
    """
    print(f"# Reading from: {filename}")
    print(f"# There are {len(mer_list)} mers(kmers/palindromes) being considered")
    print(f"# About {len_sequence - max_k + 1} positions where the mers(kmers/palindromes) could be found")
    print(f"# Reporting mers(kmers/palindromes) from length {min_k} to {max_k}")
    print(f"# which are under- or over-represented with E_value < {max_e}\n")
    print(
        f"#{'kmer'}\t{'Observed':>10}\t{'Expected':>10}\t{'Z_score':>10}\t{'Evalues':>10}\t{'Odds':>10}\t{'Diff':>10}\t{'Scores':>10}\t{'NScores':>10}\t{'Log_odds':>10}\n")
    sort_data = copy.deepcopy(data_kmer)
    sort_data = sorted(sort_data, key=operator.itemgetter(3))
    for k_data in sort_data:
        kmer, obs, exp, z_scr, e_val, odds, diff, scr, nscr, lod = k_data
        if e_val <= max_e:
            print(f"{kmer}\t{obs}\t{exp}\t{z_scr}\t{e_val}\t{odds}\t{diff}\t{scr}\t{nscr}\t{lod}\n")


def get_kmer_count_slide_window(sequence, alphabet, window, step, min_k, max_k):
    slide_mer_count = defaultdict(Counter)
    for chunk, s, e in get_chunks(sequence, window, step):
        pos = '_'.join((str(s), str(e)))
        slide_mer_count[pos].update(count_kmers(chunk, alphabet, min_k, max_k))
    return slide_mer_count


if __name__ == '__main__':
    alphabet = alphabet.iupac_dna
    max_e = 0.01
    filename = "Data/test/H_influenzae.fna.gz"
    for name, seq in fasta_parser.parse_fasta(filename):
        seq_len = len(seq) - count_umbiguous_bases(seq)
        pal_list = list(get_palindromes(alphabet, 4, 6))
        pal_counts = count_kmers(seq, alphabet, 4 - 2, 6)
        rev_strand_cnt = dict((get_reverse_complement(kmer), cnt) for kmer, cnt in pal_counts.items())
        for kmer, cnt in rev_strand_cnt.items():
            pal_counts[kmer] += cnt
        bases_freqs = base_stats(seq, alphabet, False, True)
        pal_freqs = kmers_frequencies(pal_counts)
        pal_expected = get_expected_values(pal_list,
                                           pal_counts)
        pal_zscores = get_z_scores(pal_list,
                                   pal_counts,
                                   pal_expected,
                                   seq_len)
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
                                        pal_counts)
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
        print_results_stats(filename,
                            pal_list,
                            seq_len,
                            4,
                            6,
                            max_e,
                            pal_data)
