#!usr/bin/env python
# -*- coding: utf-8 -*-


import math
from collections import defaultdict
import pandas as pd

def get_expected_values(kmer_list, kmer_counts):
    """Calculates the expected value of all palindrome kmers in a sequence.
    Inputs:
    pal_lst = list of palindromes (str)
    counts = dictionary of kmer counts
    Output:
    expected - dictionary with expected values for all palindromes.
    The expected values are calculated as:
       E(C(W)) = C(W1-Wn-1) * C(W2-Wn) / C(W2-Wn-1)
    """
    expected = defaultdict(float)
    for kmer in kmer_list:
        suf = kmer_counts[kmer[1:]]
        pre = kmer_counts[kmer[:-1]]
        mid = kmer_counts[kmer[1:-1]]
        # to catch the divide by zero error
        if mid == 0:
            expected[kmer] = 0.0
        else:
            ex = (suf * pre) / mid
            expected[kmer] = expected.get(kmer, 0.0) + ex
    return expected


def get_z_scores(kmer_list, kmer_counts, expected_kmers, len_seq):
    """Calculates the z_score of all palindromes.
    Input:
    palindrome_lst = list of palindromes (str)
    counts = dictionary of kmer counts
    expected = dictionary of kmer expected values
    length_sequence = length of sequence (int)
    Output:
    z_score dictionary where key are palindromes and values are the calculated z_score (float)
    The z_scores are calculated as:
        Z(W) = (C(W)) - E(C(W)) / sigma(W)
    And sigma as:
        sigma(W) = sqrt(E(C(W))) * (1 - E(C(W)/N))
    """
    z_score = defaultdict(float)
    for kmer in kmer_list:
        if expected_kmers[kmer] == 0.0:
            z_score[kmer] = 0.0
        else:
            sigma = math.sqrt(expected_kmers[kmer]) * (1 - expected_kmers[kmer] / (2 * len_seq))
            z = (kmer_counts[kmer] - expected_kmers[kmer]) / sigma
            z_score[kmer] = z_score.get(kmer, 0.0) + z
    return z_score


def get_pvalues(kmer_list, z_score_kmers):
    """Calculates the p_value of all palindromes.
    Input:
    palindrome_lst - list of palindromes (str)
    z_score - a dictionary with palindromes z_scores (float)
    Output:
    a palindrome p_value dictionary.
    For probability of getting a z-value larger than t
    P(z > t) = erfc(t/sqrt(2))/2
    For probability of getting a z-value smaller than t
    P(z > t) = erfc(-t/sqrt(2))/2
    """
    p_values = defaultdict(float)
    for kmer in kmer_list:
        if z_score_kmers[kmer] < 0.0:
            under = math.erfc(-z_score_kmers[kmer] / math.sqrt(2)) / 2
            p_values[kmer] = p_values.get(kmer, 0.0) + under
        else:
            over = math.erfc(z_score_kmers[kmer] / math.sqrt(2)) / 2
            p_values[kmer] = p_values.get(kmer, 0.0) + over
    return p_values


def get_evalues(kmer_list, p_value_kmers):
    """Calculates the e_value of all palindrome kmers.
    Inputs:
    palindrome_lst - list of palindromes (str)
    p_value - a dictionary where key are the palindromes (str) and the value are p_value (float)
    Output:
    The e_value as a dictionary where the key are the palindrome (str)
    and value are e_value (float).
    """
    e_value = defaultdict(float)
    num_tests = len(kmer_list)
    for kmer in kmer_list:
        p = p_value_kmers[kmer] * num_tests
        e_value[kmer] = e_value.get(kmer, 0.0) + p
    return e_value


def get_scores(kmer_list, kmer_counts, expected_kmers):
    """Calculates de kmer escore and its expected frequency in the genome.
    Inputs:
    kmer_list - list substring of length k
    kmer_counts - a dictionary with kmer and it counts
    expected_kmers - a dictionary with kmer and it expected counts
    Output:
    scores - a dictionary with kmers and it scores, calculated as:
    S = obs - exp / obs + exp
    """
    scores = defaultdict(float)
    for kmer in kmer_list:
        if expected_kmers[kmer] == 0.0 or kmer_counts[kmer] == 0.0:
            scores[kmer] = 0.0
        else:
            scr = (kmer_counts[kmer] - expected_kmers[kmer]) / (kmer_counts[kmer] + expected_kmers[kmer])
            scores[kmer] = scores.get(kmer, 0.0) + scr
    return scores


def get_new_scores(kmer_list, kmer_counts, expected_kmers):
    """Calculates de kmer escore and its expected frequency in the genome.
    Inputs:
    kmer_list - list substring of length k
    kmer_counts - a dictionary with kmer and it counts
    expected_kmers - a dictionary with kmer and it expected counts
    Output:
    scores - a dictionary with kmers and it scores, calculated as:
    S = obs - exp / obs + exp
    """
    scores = defaultdict(float)
    for kmer in kmer_list:
        scores[kmer] = scores.get(kmer, 0.0)
        if expected_kmers[kmer] == 0.0 or kmer_counts[kmer] == 0.0:
            scores[kmer] = 0.0
        else:
            scr = kmer_counts[kmer] / (kmer_counts[kmer] + expected_kmers[kmer])
            scores[kmer] = scores.get(kmer, 0.0) + scr
    return scores


def get_odds_ratio(kmer_list, kmer_freqs):
    ors = defaultdict(float)
    for kmer in kmer_list:
        midf = kmer_freqs[kmer[1:-1]]
        pref = kmer_freqs[kmer[:-1]]
        suff = kmer_freqs[kmer[1:]]
        kmf = kmer_freqs[kmer]
        # to catch the divide by zero error
        if midf == 0.0 or kmf == 0.0 or pref == 0.0 or suff == 0.0:
            ors[kmer] = ors.get(kmer, 0.0)
        else:
            od = (kmf * midf) / (pref * suff)
            ors[kmer] = ors.get(kmer, 0.0) + od
    return ors


def get_difference(kmer_list, kmer_counts, expected_kmers):
    diff = defaultdict(float)
    for kmer in kmer_list:
        d = kmer_counts[kmer] - expected_kmers[kmer]
        diff[kmer] = diff.get(kmer, 0.0) + d
    return diff


def get_log_odds(kmer_list, kmer_counts, expected_kmers):
    log_ods = defaultdict(float)
    for kmer in kmer_list:
        log_ods[kmer] = log_ods.get(kmer, 0.0)
        if kmer_counts[kmer] == 0.0 or expected_kmers[kmer] == 0.0:
            log_ods[kmer] = 0.0
        else:
            lod = math.log(kmer_counts[kmer] / expected_kmers[kmer])
            log_ods[kmer] += lod
    return log_ods


def get_kmer_statistics(kmer_list,
                        kmer_counts,
                        kmer_expected,
                        kmer_z_scores,
                        kmer_e_values,
                        kmer_odds_ratios,
                        kmer_diffs,
                        kmer_scores,
                        kmer_new_scores,
                        kmer_log_odds):
    """"""
    stats = []
    for kmer in kmer_list:
        obs = kmer_counts[kmer]
        exp = kmer_expected[kmer]
        z_scr = kmer_z_scores[kmer]
        e_val = kmer_e_values[kmer]
        odds = kmer_odds_ratios[kmer]
        diff = kmer_diffs[kmer]
        scr = kmer_scores[kmer]
        nscr = kmer_new_scores[kmer]
        lod = kmer_log_odds[kmer]
        stats.append((kmer,
                      obs,
                      exp,
                      z_scr,
                      e_val,
                      odds,
                      diff,
                      scr,
                      nscr,
                      lod))
    return stats

def get_dataframe_from_kmer_alldata(dir_out, filename, k, kmerdata):
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
    df.to_csv(f"{dir_out}/{filename}_{k}_.csv")
