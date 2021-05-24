#!usr/bin/env python
# -*- coding: utf-8 -*-

import random
import operator
from functools import reduce
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
from fasta_parser import parse_fasta, str_punctuation_strip


def prob_seq(seq, p_gc):
    # given a GC content, what is the probability
    # of getting the particular sequence

    assert (0 <= p_gc <= 1)
    # the probability of obtaining sequence seq
    # given a background gc probability of .5
    ps = []
    for char in seq:
        if char in 'CG':
            ps.append(p_gc / 2)
        elif char in 'AT':
            ps.append((1 - p_gc) / 2)
        else:
            raise ("Unexpected char: ", char)
    return reduce(operator.mul, ps, 1)


def get_shuffled_genome(filename):
    """Write a shuffled version of the sequence to a FASTA file."""
    for name, sequence in parse_fasta(filename):
        name = name.replace(',', '').split(' ')
        name = '_'.join(name[:4])
        seq = list(sequence.upper())
        # random.shuffle(seq)
        return name, ''.join(seq)


def translate_chars(s, chars, targets):
    # make a dictionary and use it
    # ex: translate_chars('mississipi',['i','s'],['e','t']) -> mettettepe
    # ex: translate_chars('mississipi','is','et') -> mettettepe
    # note only char->char substitutions are possible. For multiple
    # substring -> substring substitutions, one can use replace
    return str.translate(s, str.maketrans(chars, targets))


def get_gc_content(sequence):
    """Returns the gc content of a genome."""
    seq_len, sequence = len(sequence), sequence.upper()
    c = sequence.count('C')
    g = sequence.count('G')
    return round((c + g) / seq_len, 4)


def get_at_content(gc):
    """Returns the at content of a genome."""
    return 1 - gc


def get_at_gc_ratio(at, gc):
    """Returns the at/gc ratio of a genome."""
    return at / gc


def count_all_bases(sequence):
    bases = set(sequence)
    all_bases = defaultdict(int)
    for base in bases:
        all_bases[base] = sequence.count(base)
    return all_bases


def gc_content_sequence_window(sequence, as_overlapp=False, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. In
    overlapp windows of lenght k."""
    sequence = sequence.upper()
    seq_len = len(sequence)
    res = []
    non_overlapp = range(0, len(sequence) - k + 1, k)
    overlapp = range(0, seq_len - k + 1)
    if as_overlapp:
        for i in overlapp:
            subseq = sequence[i:i + k]
            res.append(round(((subseq.count('C') + subseq.count('G')) / len(subseq)), 4) * 100)
    else:
        for j in non_overlapp:
            subseq = sequence[j:j + k]
            res.append(round(((subseq.count('C') + subseq.count('G')) / len(subseq)), 4) * 100)
    return res


def get_sequence_from_scaffolds(filename):
    return ''.join([seq.upper() for name, seq in parse_fasta(filename)])


def codon_frequency(seq, codon_table):
    """
    Return dict of codon frequency.
    codon_table is a list of trinucleotides
    """
    empty = Counter(dict([(c, 0) for c in codon_table]))
    tmp = [seq.upper()[i:i + 3] for
           i in range(0, len(seq), 3)]
    print(tmp)
    tmp = filter(lambda x: len(x) == 3, tmp)
    return empty + Counter(tmp)


def gc_var(sequence):
    """Returns the gc variance according to a k window squence"""
    gc = get_gc_content(sequence) * 100
    gc_i = np.array(gc_content_sequence_window(sequence, as_overlapp=False, k=20))
    len_gc_i = np.shape(gc_i)[0]
    dif = gc_i - gc
    return np.log((1 / len_gc_i) * sum(abs(dif)))


def count_umbiguous_bases(sequence):
    sequence = sequence.upper()
    amb = ['N', 'R', 'Y', 'W', 'S', 'K', 'M']
    return sum({base: sequence.count(base) for base in amb}.values())


def base_stats(sequence, alphabet, as_count=False, as_dict=False):
    """Calculates de frequency or the number of bases in a sequence.
    Inputs:
    sequence - string representing the sequence
    alphabet - a alphabet (strings characters) that compound the string sequence
    as_count - boolean set as False
    as_dict - boolean set as False
    Output:
    counts - as default returns a numpy array as frequencies (floats)
    > baseFreqs(seq, 'ACGT', asCounts = False, asDict = False)
    array([0.25, 0.25, 0.25, 0.25])

    as_count - True, returns a numpy array of counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = False)
    array([2, 2, 2, 2])

    as_dict - True and as_count as default (False) returns a dictionary as bases frequencies (float)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = False, asDict = True)
    {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    as_count True and as_dict True, returns a dictionary as base counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = True)
    {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    """
    seq = sequence.upper()
    counts = np.array([seq.count(i) for i in alphabet])
    if as_count:
        freqs = counts
    else:
        freqs = counts / sum(counts * 1.0)
    if as_dict:
        return dict(zip(alphabet, freqs))
    else:
        return freqs


def get_strand_complement(sequence):
    """Returns the complement strand of the genome."""
    seq = sequence.upper()
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """Returns the reverse complement strand of the genome."""
    seq = sequence.upper()
    return get_strand_complement(seq)[::-1]


def get_chunks(sequence, window, step=1):
    """Returns a chunk of length of window_size and the end of the window size"""
    sequence = sequence.upper()
    k = len(sequence)
    for i in range(0, k - window + 1, step):
        end = i + window
        chunk = sequence[i:end]
        assert len(chunk) == window
        yield chunk, i, end


def get_sequence_skew(sequence):
    """Returns the difference between the total number of
    occurrences of G and the total number of occurrences of C in
    the first i elements of the sequence. """
    sequence = sequence.upper()
    skew = [0]
    for idx, element in enumerate(sequence):
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            skew.append(skew[idx])
    return skew


def get_minimum_skew(sequence):
    """Returns a position in a sequence minimizing the skew."""
    min_skew = []
    skew = get_sequence_skew(sequence)
    m_skew = min(skew)
    for idx in range(len(sequence) + 1):
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


def plot_base_frequency_genome(x_data, y_data, x_label, y_label):
    base_markers = {"A": "b-",
                    "C": "r-",
                    "G": "g-",
                    "T": "y-",
                    "N": "k-"}
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111)
    y_names = []
    for y in y_data:
        y_names.append(y)
        ax.plot(x_data, y_data[y], base_markers[y], label=y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    ax.legend(y_names)
    plt.grid(True)


def base_content_slide_window(sequence, name, alphabet, window, step, plot=False):
    sequence = sequence.upper()
    bases = set(alphabet.upper())
    base_freqs = defaultdict(list)
    sizes = []
    for base in bases:
        base_freqs[base] = base_freqs.get(base, [])
    for i in range(0, len(sequence) - window + 1, step):
        subseq = sequence[i:i + window]
        assert (len(subseq) == window), 'The lenght of the subsequence must have the same size of the window'
        for base in bases:
            freq = subseq.count(base) / len(subseq) * 100
            base_freqs[base].append(round(freq, 4))
        sizes.append((i + window))
    if plot:
        plot_base_frequency_genome(sizes, base_freqs, 'Genome (kb)', 'Frequencies')
        plt.title(f"Base Distribuition in {name} genome")
        plt.savefig(f"{name}_base_freq_slidewindow_plot.pdf")
    return base_freqs, sizes


def strand_stats(sequence, alphabet, start):
    alphabet = alphabet.upper()
    seq_len, seq = len(sequence), sequence.upper()
    half_gen = (seq_len // 2)
    ter = (start + half_gen)
    strand_stat = defaultdict(tuple)
    # for circular genomes
    if ter > seq_len:
        ter = ter - (seq_len + 1)
    for base in alphabet:
        base_total = seq.count(base)
        if ter > start:
            for_strand = seq[start:ter].count(base)
            rev_strand = (base_total - for_strand)
        else:
            rev_strand = seq[ter:start].count(base)
            for_strand = (base_total - rev_strand)
        dif = (for_strand - rev_strand)
        strand_stat[base] = (base_total, for_strand, rev_strand, dif)
    return strand_stat


def print_strand_stats(str_sta):
    print('   Total\tFor\tRev\tDif')
    for base, count in str_sta.items():
        print(f'{base}: {str(count[0])}\t{str(count[1])}\t{str(count[2])}\t{str(count[3])}')


def get_name_fasta(name):
    name = name.replace(',', ' ').split(' ')
    name = '_'.join(name[:4])
    return str_punctuation_strip(name)


def weighted_choice(seq):
    return random.choice(sum(([v] * wt for v, wt in seq), []))


def difreq(seq):
    """Return the count of the dinucleotides as a dictionary.
    seq = "ababcdcdabdcabvababab"
    difreq(seq)
    {'a': {'b': 7}, ab =7
    'b': {'a': 3, 'c': 1, 'd': 1, 'v': 1},
    'c': {'d': 2, 'a': 1},
    'd': {'c': 2, 'a': 1},dc =2, da =1
    'v': {'a': 1}} va = 1
    Ex. ab
    """
    counts = defaultdict(lambda: defaultdict(int))
    for a, b in zip(seq, seq[1:]):
        counts[a][b] += 1
    return dict((k, dict(v)) for k, v in counts.items())


def shuffle_difreq(seq):
    freqs = difreq(seq)

    # get the first base by the total frequency across the sequence
    shuff_seq = [weighted_choice(Counter(seq).items())]

    for i in range(1, len(seq)):
        # each following base is based of the frequency of the previous base
        # and their co-occurence in the original sequence.
        shuff_seq.append(weighted_choice(freqs[shuff_seq[-1]].items()))

    return "".join(shuff_seq)


def count_mers(sequence, alphabet, kmin, kmax):
    """Returns a dictionary like object with kmers as keys and counts
    as values.
    Inputs:
        sequence - a sequence as a string
        alphabet - a sequence as string with all 
                   the allowed charcaters in the kmer
        kmin - integer represent the minimum kmerlength 
        kmax - integer represent the maximun kmer length 
    kmer = s substring of the sequence with length kmin<k<kmax)
    """
    kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)
    kmer_from_seq = get_kmers_from_sequence(sequence, alphabet, kmin, kmax)
    counts = defaultdict(int)
    for kmer in kmer_list:
        counts[kmer] = counts.get(kmer, 0)
    for kmer in kmer_from_seq:
        if kmer not in counts:
            counts[kmer] = 1
        else:
            counts[kmer] += 1
    return counts
