__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis
from glob import glob
from itertools import product, combinations, chain, izip
import os, os.path
from Code.GeneralUtils import prots_from_path
import logging
if __name__ == '__main__':

    testfuns = ['Mutual_Info', 'OMES', 'Linkage']#, 'SBASC_McLachlan']
    fname = 'processing--%s--%s.%s.log'

    curatedfiles = glob('/hivdata/curated/MergedDir/*.aln')
    cureout = '/hivdata/curated/LinkageResults'
    largefiles = glob('/hivdata/SubtypeB/MergedDir/*.aln')
    largeout = '/hivdata/SubtypeB/LinkageResults'

    iterable = chain(combinations(curatedfiles,2),
        izip(curatedfiles, curatedfiles),
        combinations(largefiles,2),
        izip(largefiles, largefiles))

    for f1, f2 in iterable:
        for fun in testfuns:
            p1 = f1.split(os.sep)[-1].split('.')[0]
            p2 = f2.split(os.sep)[-1].split('.')[0]

            if 'curated' in f1:
                opath = cureout
            else:
                opath = largeout

            ofile = os.path.join(opath, '%s--%s.%s.res' % (p1, p2, fun))
            if not os.path.exists(ofile):
                logging.basicConfig(filename=fname % (p1, p2, fun),level=logging.DEBUG)
                logging.warning('Processing %s, %s, %s' % (f1, f2, fun))
                PredictionAnalysis(f1, f2, ofile, granular=True, limit_functions=set([fun]))
                open(ofile+'.done', 'w').write('DONE!')
