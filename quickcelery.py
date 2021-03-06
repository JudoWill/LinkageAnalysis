__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis
from glob import glob
from itertools import product, combinations, chain, izip
import os, os.path
from Code.GeneralUtils import prots_from_path
import logging
import argparse

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Distributed')
    parser.add_argument('--start-dirty', dest = 'dirty', action = 'store_true',
        default = False)
    parser.add_argument('--do-MI', dest = 'MI', action = 'store_true',
        default = False)
    parser.add_argument('--do-OMES', dest = 'OMES', action = 'store_true',
        default = False)
    parser.add_argument('--do-link', dest = 'link', action = 'store_true',
        default = False)
    parser.add_argument('--do-SBASC', dest = 'SBASC', action = 'store_true',
        default = False)
    args = parser.parse_args()

    testfuns = []
    if args.MI:
        testfuns += ['Mutual_Info']
    if args.OMES:
        testfuns += ['OMES']
    if args.link:
        testfuns += ['Linkage']
    if args.SBASC:
        testfuns += ['SBASC_McLachlan']

    c=1
    fname = 'processing--%i.log'
    while os.path.exists(fname % c):
        c+=1


    logging.basicConfig(filename=fname % c,level=logging.DEBUG)

    curatedfiles = sorted(glob('/hivdata/curated/MergedDir/*.aln'))
    cureout = '/hivdata/curated/LinkageResults'
    largefiles = sorted(glob('/hivdata/SubtypeB/MergedDir/*.aln'))
    largeout = '/hivdata/SubtypeB/LinkageResults'

    if args.dirty:
        iterable = chain(izip(largefiles, largefiles),
            combinations(largefiles,2),
            izip(curatedfiles, curatedfiles),
            combinations(curatedfiles,2))
    else:

        iterable = chain(izip(curatedfiles, curatedfiles),
            combinations(curatedfiles,2),
            izip(largefiles, largefiles),
            combinations(largefiles,2))

    for f1, f2 in iterable:
        for fun in testfuns:
            p1 = f1.split(os.sep)[-1].split('.')[0]
            p2 = f2.split(os.sep)[-1].split('.')[0]

            if 'curated' in f1:
                opath = cureout
            else:
                opath = largeout

            ofile = os.path.join(opath, '%s--%s.%s.res' % (p1, p2, fun))
            if not (os.path.exists(ofile+'.p') or os.path.exists(ofile+'.done')):
                open(ofile+'.p', 'w').write('checking!')
                logging.warning('Processing %s, %s, %s' % (f1, f2, fun))
                if os.path.exists(ofile):
                    mode = 'a'
                else:
                    mode = 'w'
                PredictionAnalysis(f1, f2, ofile, open_mode = mode, limit_functions=set([fun]))
                open(ofile+'.done', 'w').write('DONE!')
