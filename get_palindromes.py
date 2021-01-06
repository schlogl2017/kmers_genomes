#!usr/bin/env python
# -*- coding: utf-8 -*-
import itertools
import alphabet
from sequence_utils import get_reverse_complement
import fasta_parser



def get_palindromes(alphabet, min_k, max_k):
    """Generates all DNA palindromes over the range from min_k to max_k.
    Inputs:
    min_k - minimum palindrome length (int)
    max_k - maximum palindrome length (int)
    Output:
    yields all possible DNA palindromes (str) of length min_k to max_k.
    Some definitions:
    A palindrome is defined as a sequence which is equal to its reverse-complement.
    Note: for odd length palindromes, the middle base does not need to be the same
    in the reverse-complement.
    Ex.: AAT is a legal palindrome even though its reverse-complement is ATT
    """
    for k in range(min_k, (max_k + 1)):
        for mer in itertools.product(alphabet, repeat=int(k / 2)):
            kmer = ''.join(mer)
            # even pal
            if k % 2 == 0:
                pal = kmer + get_reverse_complement(kmer)
                yield pal
            else:
                for base in alphabet:  # odd pal
                    pal = kmer + base + get_reverse_complement(kmer)
                    yield pal

if __name__ == '__main__':
    alphabet = alphabet.iupac_dna
    filename = "Data/test/H_influenzae.fna.gz"
    for name, seq in fasta_parser.parse_fasta(filename):
        pal_list = list(get_palindromes(alphabet, 4, 6))
        print(len(pal_list))
        print(pal_list)
