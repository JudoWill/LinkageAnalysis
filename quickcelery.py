__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis
from glob import glob
from itertools import product
import os, os.path
from Code.GeneralUtils import prots_from_path

if __name__ == '__main__':

    if os.path.exists('/hivdata/curated/MergedDir/'):
        ifiles = glob('/hivdata/curated/MergedDir/*.aln')
        opath = '/hivdata/curated/LinkageResults'
    else:
        ifiles = glob('HIVData/curated/MergedDir/*.aln')
        opath = 'HIVData/curated/LinkageResults'
    print ifiles
    iterable = product(ifiles, repeat=2)
    iterable.next() #env-env seems to have an issue

    for f1, f2 in product(ifiles, repeat=2):
        print 'Processing', f1, f2
        p1 = f1.split(os.sep)[-1].split('.')[0]
        p2 = f2.split(os.sep)[-1].split('.')[0]

        ofile = os.path.join(opath, '%s--%s.res' % (p1, p2))
        PredictionAnalysis(f1, f2, ofile, granular=True)


    if os.path.exists('/hivdata/SubtypeB/MergedDir/'):
        ifiles = glob('/hivdata/SubtypeB/MergedDir/*.aln')
        opath = '/hivdata/SubtypeB/LinkageResults'
    else:
        ifiles = glob('HIVData/SubtypeB/MergedDir/*.aln')
        opath = 'HIVData/SubtypeB/LinkageResults'

    for f1, f2 in product(ifiles, repeat=2):
        print 'Processing', f1, f2
        p1 = f1.split(os.sep)[-1].split('.')[0]
        p2 = f2.split(os.sep)[-1].split('.')[0]

        ofile = os.path.join(opath, '%s--%s.res' % (p1, p2))
        PredictionAnalysis(f1, f2, ofile, granular=True)