__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis
from glob import glob
from itertools import product, combinations, chain, izip
import os, os.path
from Code.GeneralUtils import prots_from_path
import logging
if __name__ == '__main__':

    fname = 'processing%i.log'
    c=1
    while not os.path.exists(fname % c):
        c+=1
    logging.basicConfig(filename=fname % c,level=logging.DEBUG)

    if os.path.exists('/hivdata/curated/MergedDir/'):
        ifiles = glob('/hivdata/curated/MergedDir/*.aln')
        opath = '/hivdata/curated/LinkageResults'
    else:
        ifiles = glob('HIVData/curated/MergedDir/*.aln')
        opath = 'HIVData/curated/LinkageResults'
    print ifiles
    iterable = chain(combinations(ifiles,2), izip(ifiles, ifiles))

    for f1, f2 in iterable:
        print 'Processing', f1, f2
        p1 = f1.split(os.sep)[-1].split('.')[0]
        p2 = f2.split(os.sep)[-1].split('.')[0]

        ofile = os.path.join(opath, '%s--%s.res' % (p1, p2))
        if not os.path.exists(ofile):
            PredictionAnalysis(f1, f2, ofile, granular=True)


    if os.path.exists('/hivdata/SubtypeB/MergedDir/'):
        ifiles = glob('/hivdata/SubtypeB/MergedDir/*.aln')
        opath = '/hivdata/SubtypeB/LinkageResults'
    else:
        ifiles = glob('HIVData/SubtypeB/MergedDir/*.aln')
        opath = 'HIVData/SubtypeB/LinkageResults'

    print ifiles
    iterable = chain(combinations(ifiles,2), izip(ifiles, ifiles))

    for f1, f2 in iterable:
        print 'Processing', f1, f2
        p1 = f1.split(os.sep)[-1].split('.')[0]
        p2 = f2.split(os.sep)[-1].split('.')[0]

        ofile = os.path.join(opath, '%s--%s.res' % (p1, p2))
        if not os.path.exists(ofile):
            PredictionAnalysis(f1, f2, ofile, granular=True)