#!usr/bin/env python
# -*- coding: utf-8 -*-



def get_log_odds(kmer_list, kmer_counts, expected_kmers):
    log_ods = defaultdict(float)
    for kmer in kmer_list:
        if kmer_counts[kmer] == 0.0 or expected_kmers[kmer] == 0.0:
            log_ods[kmer] = 0.0
        else:
            lod = math.log(kmer_counts[kmer] / expected_kmers[kmer])
            log_ods[kmer] = log_ods.get(kmer, 0) + lod
    return log_ods



