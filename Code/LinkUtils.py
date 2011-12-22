from __future__ import division
__author__ = 'will'


from operators import gt, lt
from random import shuffle


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