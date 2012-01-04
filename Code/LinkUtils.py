from __future__ import division
from itertools import product
from math import sqrt, log
from operator import itemgetter
import glob, os, os.path
__author__ = 'will'

from functools import partial
from operator import gt, lt
from random import shuffle
from collections import defaultdict
from pylru import lrudecorator
from celery.task import task
from celery.exceptions import TimeoutError, WorkerLostError, SoftTimeLimitExceeded, TimeLimitExceeded
from Queue import Queue, Empty
import logging

class LinkCalculator(object):

    def __init__(self):
        sub_mats = glob.glob('Code/*.mat')
        self.process_funcs = []
        for f in sub_mats:
            name = os.path.splitext(os.path.basename(f))[0]
            self.process_funcs.append(('SBASC_'+name,
                partial(calculate_SBASC, get_sub_mat(f))))

        self.process_funcs.append(('Mutual_Info', calculate_mutual_info))
        self.process_funcs.append(('OMES', calculate_OMES))
        self.process_funcs.append(('Linkage', calculate_mapping))
        self.sufs = ['_raw', '_pval', '_null', '_count']

    def calculate_all(self, seq1, seq2):

        for bname, func in self.process_funcs:
            results = calculate_vals(list(seq1), list(seq2), func)
            for suffix, result in zip(self.sufs, results):
                yield bname+suffix, result

    def get_fields(self):

        fields = []
        for bname, _ in self.process_funcs:
            for suffix in self.sufs:
                fields.append(bname+suffix)
        return fields

def check_precision(extreme_counts, total_counts, logprecision):
    """ Returns True if there is enogh precision to accept the p-value
    """
    extreme_counts = max(extreme_counts, 1)
    total_counts = max(total_counts, 1)
    return (log(total_counts,10) + log(extreme_counts/total_counts,10)) > logprecision

def function_linker(func, ls1, ls2, preargs, distributed = True, **kwargs):
    if len(preargs)==1:
        args = (preargs[0], ls1, ls2)
    else:
        args = (ls1, ls2)
    if distributed:
        return func.delay(*args, **kwargs)
    else:
        return func(*args, **kwargs)


def calculate_vals(s1, s2, testfun, key = gt, preargs = (), minreps = 500, maxreps = 1e6):

    ls1 = list(s1)
    ls2 = list(s2)
    extreme_count = 0
    total_sum = 0
    total_count = 0

    try:
        trueval = function_linker(testfun, ls1, ls2, preargs, distributed=False)
    except ZeroDivisionError:
        return 0, 1.0, 0, 0

    reslist = function_linker(testfun, ls1, ls2, preargs, distributed=False,
        shuf = True, batch = minreps)
    for res in reslist:
        if res > trueval:
            extreme_count += 1
        total_count += 1
        total_sum += res

    while total_count < maxreps and not check_precision(extreme_count, total_count, 1):
        reslist = function_linker(testfun, ls1, ls2, preargs, distributed=False,
            shuf = True, batch = minreps)
        for res in reslist:
            if res > trueval:
                extreme_count += 1
            total_count += 1
            total_sum += res
        logging.info('pval: %f, mean: %f, count: %i' % (extreme_count/total_count, total_sum/total_count, total_count))


    return trueval, extreme_count/total_count, total_sum/total_count, total_count




def celery_calculate_vals(s1, s2, testfun, preargs = (), key = gt, minreps = 500, maxreps = 1e6):



    def process_que(que, trueval, timeout):
        qsize = que.qsize()
        total_count = 0
        true_count = 0
        total_sum = 0.0
        flag = False
        for num in xrange(qsize):
            logging.info('getting que num %i' % num)
            try:
                asyncres = que.get_nowait()
                reslist = asyncres.get(timeout=timeout)
                for res in reslist:
                    if res > trueval:
                        true_count += 1
                    total_count += 1
                    total_sum += res
            except Empty:
                #queue has nothing left
                break
            except TimeoutError:
                logging.warning('putting back in')
                que.put(asyncres)
            except WorkerLostError:
                pass
            except TimeLimitExceeded:
                flag = True
            except:
                logging.error('something else!')
        return true_count, total_sum, total_count, flag

    batchsize = 1000
    groupingsize = 20
    ibatch = int((minreps*1.5)/groupingsize)+1 #enough so if even HALF don't finish it will still pass minreps
    ls1 = list(s1)
    ls2 = list(s2)
    que = Queue()

    try:
        trueval = function_linker(testfun, ls1, ls2, preargs, distributed=False)
    except ZeroDivisionError:
        return 0, 1.0, 0, 0

    logging.warning('doing initial batch')
    for _ in xrange(groupingsize):
        que.put(function_linker(testfun, ls1, ls2, preargs, shuf = True, batch = ibatch))
    extreme_count, total_sum, total_count, flag = process_que(que, trueval, 0.1)

    if total_count >= minreps and check_precision(extreme_count, total_count, 1):
        logging.info('Initial batch was enough!')
    else:
        while total_count < maxreps and not check_precision(extreme_count, total_count, 1):
            logging.info('putting in %i' % int(batchsize))
            for _ in xrange(groupingsize):
                que.put(function_linker(testfun, ls1, ls2, preargs, shuf = True, batch = int(batchsize)))
            e, s, t, flag = process_que(que, trueval, 5)
            if flag:
                batchsize = max(int(batchsize*0.5), ibatch)
            extreme_count += e
            total_sum += s
            total_count += t

            logging.info('pval: %f, mean: %f, count: %i' % (extreme_count/total_count, total_sum/total_count, total_count))

    logging.info('emptying que')
    e, s, t, _ = process_que(que, trueval, 0.1)
    extreme_count += e
    total_sum += s
    total_count += t

    return extreme_count/total_count, total_sum/total_count, total_count


def signal2prob(signal):
    counts = signal2count(signal)
    num = len(signal)
    for key, val in counts.items():
        counts[key] = val/num
    return counts

def signal2count(signal):
    counts = defaultdict(int)
    for s in signal:
        counts[s] += 1
    return counts



@task()
def calculate_mutual_info(signal1, signal2, shuf = False, batch=False, **kwargs):
    """Caluculates the Mutual Information shared by two signals.

    Arguements:
    signal1 -- An iterable indicating the first signal
    signal2 -- An iterable indicating the second signal

    Signals MUST be the same length! Items must be hashable!

    Returns:
    Mutual Information -- float"""



    if shuf:
        shuffle(signal1)
        shuffle(signal2)

    if batch:
        res = []
        extra = {}
        try:
            for _ in xrange(batch):
                r, extra = calculate_mutual_info(signal1, signal2, shuf=True, want_extra=True, **extra)
                res.append(r)
        except SoftTimeLimitExceeded:
            pass
        return res

    overlap_prob = signal2prob(zip(signal1, signal2))
    signal1_prob = kwargs.get('S1prob',signal2prob(signal1))
    signal2_prob = kwargs.get('S2prob',signal2prob(signal2))

    num_items = len(signal1)
    mut_info = float()


    for (s1, s2), count in overlap_prob.items():
        mut_info += overlap_prob[(s1, s2)]*log(overlap_prob[(s1, s2)]/(signal1_prob[s1]*signal2_prob[s2]))

    if kwargs.get('want_extra', False):
        return mut_info, {'S1prob':signal1_prob, 'S2prob':signal2_prob}
    else:
        return mut_info

@task()
def calculate_PNAS(signal1, signal2, shuf = False, batch = False):

    if shuf:
        shuffle(signal1)
        shuffle(signal2)


    if batch:
        res = []
        try:
            for _ in xrange(batch):
                res.append(calculate_PNAS(signal1, signal2, shuf=True))
        except SoftTimeLimitExceeded:
            pass
        return res

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

@task()
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


    counts = signal2count(zip(signal1, signal2))

    mapping = []
    while counts:
        (s1, s2), val = max(counts.items(), key = itemgetter(1))
        mapping.append((s1,s2,val))
        for ks1, ks2 in counts.keys():
            if ks1 == s1:
                counts.pop((ks1, ks2))
    return mapping

@task()
def calculate_mapping(signal1, signal2, shuf = False, batch = False):


    if batch:
        res = []
        try:
            for _ in xrange(batch):
                res.append(calculate_mapping(signal1, signal2, shuf=True))
        except SoftTimeLimitExceeded:
            pass
        return res

    if shuf:
        shuffle(signal1)
        shuffle(signal2)

    res = prediction_mapping(signal1, signal2)
    return sum(r[2] for r in res)/len(signal1)

@task()
def calculate_OMES(signal1, signal2, shuf = False, batch = False, **kwargs):
    """Finds the Observed Minus Expected Squared score.

    Calcualtes the score using the formula:
    sum((Nobs-Nex)^2/Nvalid) for all pairs in S1,S2

    Reference: Fodor et. al. 2004, PMID: 15211506
    """

    if shuf:
        shuffle(signal1)
        shuffle(signal2)

    if batch:
        res = []
        extra = {}
        try:
            for _ in xrange(batch):
                r, extra = calculate_OMES(signal1, signal2, shuf=True, want_extra=True, **extra)
                res.append(r)
        except SoftTimeLimitExceeded:
            pass
        return res

    s1_counts = kwargs.get('S1counts', signal2count(signal1))
    s2_counts = kwargs.get('S2counts', signal2count(signal2))
    Nobs = signal2count(zip(signal1, signal2))
    Nvalid = 0

    omes = 0.0
    for (s1, s2), obs in Nobs.items():
        Nex = s1_counts[s1]*s2_counts[s2]/Nvalid
        omes += ((obs-Nex)**2)/Nvalid

    if kwargs.get('want_extra', False):
        return omes, {'S1counts': s1_counts, 'S2counts':s2_counts}
    else:
        return omes

@task()
def calculate_SBASC(sub_mat, signal1, signal2, shuf = False, batch = False, **kwargs):
    """Calculates the Substitutions Based Correlation

     Determines the correlation of mutations based on any substituion matrix.

     Based on Olmea et. al. 1999 ... PMID: 10547297

    """

    if shuf:
        shuffle(signal1)
        shuffle(signal2)

    if batch:
        res = []
        extra = {}
        try:
            for _ in xrange(batch):
                r, extra = calculate_SBASC(sub_mat, signal1, signal2, shuf=True, want_extra=True, **extra)
                res.append(r)
        except SoftTimeLimitExceeded:
            pass
        return res

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

    s1scores, s1std = kwargs.get('S1sub', sub_score(signal1, sub_mat))
    s2scores, s2std = kwargs.get('S2sub', sub_score(signal2, sub_mat))

    #print s1scores, s1std, s2scores, s2std
    denom = s1std*s2std
    McBASC = 0.0

    for s1, s2 in zip(s1scores, s2scores):
        McBASC += s1*s2/denom
    McBASC /= N2

    if kwargs.get('want_extra', False):
        odict = {'S1sub': (s1scores, s1std),
                 'S2sub': (s2scores, s2std)}
        return McBASC, odict
    else:

        return McBASC

def get_all_sub_mats():

    sub_mats = glob.glob('Code/*.mat')
    output = []
    for f in sub_mats:
        name = os.path.splitext(os.path.basename(f))[0]
        mat = get_sub_mat(f)
        output.append((name,mat))
    return output


def get_sub_mat(filename):

    submat = defaultdict(int)
    with open(filename) as handle:
        for line in handle:
            parts = line.strip().split('\t')
            submat[(parts[0], parts[1])] = float(parts[2])
    return submat