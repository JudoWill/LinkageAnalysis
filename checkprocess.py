__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis
import os, os.path


if __name__ == '__main__':

    basedir = '/hivdata/test/'
    testfuns = ['Mutual_Info', 'OMES', 'Linkage', 'SBASC_McLachlan']
    basealign = '/hivdata/test/Env.aln'
    for tfun in testfuns:
        ofile = '/hivdata/test/Env--Env.aln.'+tfun
        if not os.path.exists(ofile):
            PredictionAnalysis(basealign, basealign, ofile, granular=True, limit_functions=set([tfun]))
