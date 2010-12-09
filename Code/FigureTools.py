from __future__ import division
import matplotlib
matplotlib.use('Agg')

import csv, os.path, os
from itertools import product, groupby
from matplotlib import axes, pyplot
from copy import deepcopy
from operator import eq, lt, gt, itemgetter
import numpy

def guessing_figures(load_dir):
    
    min_rows = 20
    rows = []
    proteins = set()
    structures = set()
    for f in os.listdir(load_dir):
        if f.endswith('.res'):
            parts = f.split('.')[0].split('--')
            proteins.add(parts[0])
            structures.add(parts[1])
            path = os.path.join(load_dir, f)
            with open(path) as handle:
                rows += list(csv.DictReader(handle, delimiter = '\t'))
    for row in rows:
        row['Linkage'] = float(row['Linkage'])
        row['3dDist'] = float(row['3dDist'])
    
    mdict = dict()
    for key, srows in groupby(rows, itemgetter('Structure')):
        mdict[key] = max([x['3dDist'] for x in srows])
        
    for row in rows:
        row['3dNorm'] = float(row['3dDist']/mdict[row['Structure']])
        
    dist_cuts = [None] + range(2,60,2)
    link_cuts = [None] + [x/10 for x in range(11)]
    proteins.add(None)
    structures.add(None)
    order = ('3dDist', 'Linkage', 'Protein', 'Structure', 'Normed')
    funs = (lt, gt, eq, eq)
    fields = order + ('FigNum','Corr')
    handle = csv.DictWriter(open(os.path.join(load_dir, 'figures', 'key.txt'), 'w'),
                            fieldnames = fields, delimiter = '\t')
    count = 1
    for args in product(dist_cuts, link_cuts, proteins, structures, [True, False]):
        trows = deepcopy(rows)
        for arg, key, test in zip(args[:-1], order, funs):
            if arg is not None:
                trows = [x for x in trows if test(x[key], arg)]
        if len(trows) > min_rows:
            
            odict = dict(zip(order, args))
            odict['FigNum'] = count
            odict['Normed'] = args[-1]
            
            if odict['Normed']:
                X = [x['3dNorm'] for x in trows]
            else:
                X = [x['3dDist'] for x in trows]
            Y = [x['Linkage'] for x in trows]
            cor = numpy.corrcoef(numpy.array([X,Y]))[1,0]
            odict['Corr'] = cor
            if cor > 0.7:
                fig, ax = pyplot.subplots(1)
                ax.scatter(X, Y)
                fig.savefig(os.path.join(load_dir, 'figures','Figure%i.png' % count))
                del(fig)
            
                print odict.items()
                handle.writerow(odict)
                count += 1
            else:
                print 'Not', odict.items()