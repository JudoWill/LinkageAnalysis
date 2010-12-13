from __future__ import division
import matplotlib
matplotlib.use('Agg')

import csv, os.path, os, tempfile
from itertools import product, groupby
from matplotlib import axes, pyplot
from copy import deepcopy
from operator import eq, lt, gt, itemgetter
import numpy



class CircosGraphs():
    
    def __init__(self):
        self.genome = list()
        self.links = list()

    def add_links(self, iterable):
        
        if self.links is None:
            self.links = tuple(iterable)
        else:
            self.links += tuple(iterable)

    def add_gene(self, name, length, color = None):
        self.genome.append({'name':name, 'color':color, 
                            'length':length})

    def genes_from_file(self, fname = 'Code/genome_template.tmp'):
        genomes = 
        with open(fname) as handle:
            for line in handle:
                if line.startswith('chr - '):
                    parts = line.strip().split()
                    self.add_gene(parts[1]. parts[4], color = parts[-1])
                    


    def make_figure(self, image_path, info_dict, link_filter = None, SCRATH = '/tmp/'):
        
        image_path = os.path.abspath(image_path)
        image_direc = image_path.rsplit(os.sep)[0]
        image_name = image_path.rsplit(os.sep)[1]
        defaults = {'kayrotype':'kayrotype.txt',
                    'outdir':image_direc,
                    'outfile':image_name,
                    'infile':'links.txt',
                    'radius':5000,
                    }

        circos_template = 'Code/circos.tem'
        source_ktype = 'Code/hiv_genome.txt'
        tmp_loc = tempfile.mkdtemp(dir = SCRATCH)
        link_file = os.path.join(tmp_loc, 'links.txt')
        circos_file = os.path.join(tmp_loc, 'circos.conf')
        ktype_file = os.path.join(tmp_loc, 'kayrotype.txt')

        with open(link_file, 'w') as handle:
            for i, link in enumerate(self.links):
                tdict = dict(**link)
                tdict['lnum'] = i
                handle.write('hivlink%(lnum)i %(chrom)s %(start)i %(stop)i\n' % tdict)
                handle.write('hivlink%(lnum)i %(chrom)s %(start)i %(stop)i\n' % tdict)

        with open(circos_template) as handle:
            ctemp = handle.read()


        with open(circos_file, 'w') as handle:
            handle.write(ctemp % defaults)
        
        shutil.copy(source_ktype, ktype_file)
        
        
        
        




def guessing_figures(load_dir):
    
    min_rows = 50
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
        row['dDistL'] = row['3dDist']
        row['dDistU'] = row['3dDist']
    
    mdict = dict()
    for key, srows in groupby(rows, itemgetter('Structure')):
        mdict[key] = max([x['3dDist'] for x in srows])
        
    for row in rows:
        row['3dNorm'] = float(row['3dDist']/mdict[row['Structure']])
        
    up_cuts = [None] + range(2,60,2)
    low_cuts = [None] + range(2,60,2)
    link_cuts = [None] + [x/10 for x in range(11)]
    proteins.add(None)
    structures.add(None)
    order = ('dDistU', 'dDistL','Linkage', 'Protein', 'Structure', 'Normed')
    funs = (lt, gt, gt, eq, eq)
    fields = order + ('FigNum','Corr')
    handle = csv.DictWriter(open(os.path.join(load_dir, 'figures', 'key.txt'), 'w'),
                            fieldnames = fields, delimiter = '\t')
    ghandle = csv.DictWriter(open(os.path.join(load_dir, 'figures', 'gkey.txt'), 'w'),
                            fieldnames = fields, delimiter = '\t')
    count = 1
    for args in product(up_cuts, low_cuts, link_cuts, proteins, structures, [True, False]):
        trows = rows
        for arg, key, test in zip(args[:-1], order, funs):
            if arg is not None:
                trows = [x for x in trows if test(x[key], arg)]
        if len(trows) > min_rows:
            
            odict = dict(zip(order, args))
            odict['FigNum'] = count
            odict['Normed'] = args[-1]
            for key, val in odict.items():
                if val is None:
                    odict[key] = 'All'
            
            if odict['Normed']:
                X = [x['3dNorm'] for x in trows]
            else:
                X = [x['3dDist'] for x in trows]
            Y = [x['Linkage'] for x in trows]
            cor = numpy.corrcoef(numpy.array([X,Y]))[1,0]
            odict['Corr'] = cor
            if abs(cor) > 0.65:
                fig, ax = pyplot.subplots(1)
                ax.scatter(X, Y)
                for key, val in odict.items():
                    odict[key] = str(val)
                
                tstr = '%(dDistL)s < 3d '
                tstr += '< %(dDistU)s ; '
                tstr += 'Link > %(Linkage)s ; '
                tstr += 'Protien: %(Protein)s ; '
                tstr += 'Structure: %(Structure)s ; '
                tstr += 'Cor: %(Corr)s'
                tstr = tstr % odict
                ax.set_title(tstr)
                fig.savefig(os.path.join(load_dir, 'figures','Figure%i.png' % count))
                del(fig)
                ghandle.writerow(odict)
                print tstr
                count += 1
            else:
                odict['FigNum'] = 0
                
            handle.writerow(odict)
