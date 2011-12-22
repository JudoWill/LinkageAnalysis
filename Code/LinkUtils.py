from __future__ import division
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