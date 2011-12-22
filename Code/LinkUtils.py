from __future__ import division
from itertools import product
from math import sqrt
from operator import itemgetter

__author__ = 'will'


from operators import gt, lt
from random import shuffle
from collections import defaultdict


def calculate_vals(s1, s2, testfun, key = gt, minreps = 500, maxreps = 1e6):

    tcount = 0
    count = 0
    total = 0
    mcut = 0.01
    trueval = testfun(s1, s2)

    while tcount < minreps or (tcount < maxreps and count/tcount > mcut):
        tcount += 1
        res = testfun(shuffle(s1), shuffle(s2))
        total += res
        if key(res, trueval):
            count += 1

    return trueval, count/tcount, total/tcount, tcount


def calculate_mutual_info(signal1, signal2, **kwargs):
    """Caluculates the Mutual Information shared by two signals.

    Arguements:
    signal1 -- An iterable indicating the first signal
    signal2 -- An iterable indicating the second signal

    Signals MUST be the same length! Items must be hashable!

    Returns:
    Mutual Information -- float"""

    def count2prob(d, num):
        for key, val in d.items():
            d[key] = val/num
        return d


    overlap = defaultdict(int)
    num_items = len(signal1)
    signal1_hist = defaultdict(int)
    signal2_hist = defaultdict(int)

    for s1, s2 in zip(signal1, signal2):
        overlap[(s1, s2)] += 1
        signal1_hist[s1] += 1
        signal2_hist[s2] += 1

    mut_info = float()
    overlap_prob = count2prob(overlap, num_items)
    signal1_prob = count2prob(signal1_hist, num_items)
    signal2_prob = count2prob(signal2_hist, num_items)

    for (s1, s2), count in overlap.items():
        mut_info += overlap_prob[(s1, s2)]*log(overlap_prob[(s1, s2)]/(signal1_prob[s1]*signal2_prob[s2]))

    return mut_info


def calculate_PNAS(signal1, signal2, **kwargs):


    snum = len(signal1)
    s1_counts = defaultdict(int)
    s2_counts = defaultdict(int)

    for s1, s2 in zip(signal1, signal2):
        s1_counts[s1] += 1
        s2_counts[s2] += 1

    s1_cons, _ = max(s1_counts.items(), key = itemgetter(1))
    s2_cons, _ = max(s2_counts.items(), key = itemgetter(1))

    s1_bin = [s1 != s1_cons for s1 in signal1] #True if MUTATION!!!
    s2_bin = [s2 != s2_cons for s2 in signal2]

    f1 = sum(s1_bin)
    f2 = sum(s2_bin)
    f12 = sum(1 for s1, s2 in zip(s1_bin, s2_bin) if s1 and s2)

    f1m = f1/snum
    f2m = f2/snum

    V1 = sum((f1m-x)**2 for x in s1_bin)
    V2 = sum((f2m-x)**2 for x in s2_bin)

    return (f12/snum - f1m*f2m)/sqrt(V1*V2)


def prediction_mapping(signal1, signal2):
    """Calculates the mapping between any two signals.

    Uses a depth-first search algorithm to match the most likely value in
    signal2 for each value in signal1. Returns a nested list of the mappings
    and thier occurance values.

    Arguements:
    signal1 -- An iterable indicating the first signal
    signal2 -- An iterable indicating the second signal

    Signals MUST be the same length! Items must be hashable!

    Returns:
    [(s1a, s2a, #occurance), (s1b, s2b, #occurance), ...]"""



    counts = defaultdict(int)
    for s1, s2 in zip(signal1, signal2):
        counts[(s1, s2)] += 1

    mapping = []
    while counts:
        (s1, s2), val = max(counts.items(), key = itemgetter(1))
        mapping.append((s1,s2,val))
        for ks1, ks2 in counts.keys():
            if ks1 == s1:
                counts.pop((ks1, ks2))
    return mapping


def calculate_OMES(signal1, signal2, **kwargs):
    """Finds the Observed Minus Expected Squared score.

    Calcualtes the score using the formula:
    sum((Nobs-Nex)^2/Nvalid) for all pairs in S1,S2

    Reference: Fodor et. al. 2004, PMID: 15211506
    """

    s1_counts = defaultdict(int)
    s2_counts = defaultdict(int)
    Nobs = defaultdict(int)
    Nvalid = 0

    for s1, s2 in zip(signal1, signal2):
        if s1 != '-' and s2 != '-':
            s1_counts[s1] += 1
            s2_counts[s2] += 1
            Nobs[(s1,s2)] += 1
            Nvalid += 1

    omes = 0.0
    for (s1, s2), obs in Nobs.items():
        Nex = s1_counts[s1]*s2_counts[s2]/Nvalid
        omes += ((obs-Nex)**2)/Nvalid

    return omes


def calculate_SBASC(sub_mat, signal1, signal2, **kwargs):
    """Calculates the Substitutions Based Correlation

     Determines the correlation of mutations based on any substituion matrix.

     Based on Olmea et. al. 1999 ... PMID: 10547297

    """

    def sub_score(signal, sub_mat):
        scores = []
        for i,j in product(xrange(len(signal)), repeat=2):
            if i != j:
                scores.append(sub_mat[(signal[i], signal[j])])
        mscore = sum(scores)/len(scores)
        stdscore =  sqrt(sum((s-mscore)**2 for s in scores)/len(scores))
        cscore = [s - mscore for s in scores]
        return cscore, stdscore

    N2 = len(signal1)**2

    s1scores, s1std = sub_score(signal1, sub_mat)
    s2scores, s2std = sub_score(signal2, sub_mat)

    denom = s2std*s2std
    McBASC = 0.0

    for s1, s2 in zip(s1scores, s2scores):
        McBASC += s1*s2/denom
    McBASC = McBASC/N2

    return McBASC


def get_sub_mat(filename):

    submat = defaultdict(int)
    with open(filename) as handle:
        for line in handle:
            parts = line.strip().split('\t')
            submat[(parts[0], parts[1])] = float(parts[2])
    return submat