from __future__ import division
import time
from math import log, exp
#from mpmath import loggamma
from operator import itemgetter, ge
from functools import partial
from itertools import imap, starmap, repeat
from random import shuffle
from Code.AlignUtils import prediction_mapping

try:
    from memorised.decorators import memorise
except ImportError:
    class memorise(object):
        def __init__(self, func):
            self.func = func
            self.cache = {}
        def __call__(self, *args, **kwargs):
            try:
                return self.cache[args]
            except KeyError:
                v = self.func(*args, **kwargs)
                self.cache[args] = v
                return v
            except TypeError:
                print 'unhash-able'
                return self.func(*args, **kwargs)
        def __repr__(self):
            return self.func.__doc__
        def __get__(self, obj, objtype):
            return partial(self.__call__, obj)

def repeatfunc(func, times , *args):
    print times
    print args
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))

def tuple_shuffle(tup):
    t = list(tup)
    shuffle(t)
    return tuple(t)

@memorise
def linkage_pval(sa, sb, num_reps = 100000):

    def get_score(signal_a, signal_b):
        mappings = prediction_mapping(tuple(signal_a), tuple_shuffle(signal_b))
        return sum(imap(itemgetter(2), mappings))/len(signal_a)
    
    rscore = get_score(sa, sb)
    check_score = partial(ge, rscore)
    return sum(imap(check_score, repeatfunc(get_score, num_reps, sa, sb)))/num_reps


PVAL_CUT = 0.05
@memorise
def logchoose(ni, ki):
    #n = max(ni, ki)
    #k = min(ni, ki)
    try:
        lgn1 = loggamma(ni+1)
        lgk1 = loggamma(ki+1)
        lgnk1 = loggamma(ni-ki+1)
    except ValueError:
        #print ni,ki
        raise ValueError


    return lgn1 - (lgnk1 + lgk1)

@memorise
def gauss_hypergeom(X, n, m, N):
    """Returns the probability of drawing X successes of m marked items
     in n draws from a bin of N total items."""

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)


    r1 = logchoose(m, X)
    try:
        r2 = logchoose(N-m, n-X)
    except ValueError:
        return 0
    r3 = logchoose(N,n)

    return exp(r1 + r2 - r3)

@memorise
def hypergeo_cdf(X, n, m, N):

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

    s = 0
    for i in range(1, X+1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s,0.0), 1)
