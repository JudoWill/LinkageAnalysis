from __future__ import division
import matplotlib
matplotlib.use('Agg')

from Code.GeneralUtils import *
from Code.AlignUtils import *
from Code.StatUtils import *
import csv, os.path, os, tempfile, shutil
from itertools import product, groupby, islice, ifilter
from matplotlib import axes, pyplot
from copy import deepcopy
from operator import eq, lt, gt, itemgetter
import numpy
from subprocess import call
import shlex
from operator import itemgetter



class CircosGraph():
    
    def __init__(self):
        self.genome = list()
        self.links = list()
    
    @staticmethod
    def load_from_dir(load_dir, align_dir):
        
        def treatrow(row, sprot, tprot, smapping, tmapping):
            nrow = dict(**row)
            int_items = ('Source-Start', 'Source-End', 
                          'Target-Start', 'Target-End',
                            'Correct-Num', 'Total-Num')
            for it in int_items:
                nrow[it] = int(nrow[it])
            tmp = (('Source-Start', smapping),
                    ('Source-End', smapping),
                    ('Target-Start', tmapping),
                    ('Target-End', tmapping))
            for it, mapping in tmp:
                try:
                    nrow[it] = mapping[nrow[it]]
                except:
                    print sprot, tprot, nrow[it], len(smapping)
                    raise IndexError

            nrow['Score'] = float(nrow['Total-Score'])
            nrow['Source-Prot'] = sprot
            nrow['Target-Prot'] = tprot
            return nrow

        ref = 'K03455'
        mapping_dict = {}
        for f in os.listdir(align_dir):
            if f.endswith('.aln'):
                key = f.split('.')[0]
                align = Alignment.alignment_from_file(os.path.join(align_dir, f))

                _, mapping_dict[key] = align.convert_numbering(ref)
                print key, len(mapping_dict[key])

        grouper = itemgetter('Source-Start', 'Source-End', 
                                'Target-Start', 'Target-End')

        graph = CircosGraph()
        graph.genes_from_file()
        for f in os.listdir(load_dir):
            if f.endswith('.res'):
                parts = f.split('.')[0].split('--')
                source = parts[0]
                target = parts[1]
                smapping = mapping_dict[source]
                tmapping = mapping_dict[target]
                with open(os.path.join(load_dir, f)) as handle:
                    for key, rows in groupby(csv.DictReader(handle, delimiter = '\t'), grouper):
                        row = rows.next()
                        if row['Correct-Num'] is not None and not row['Correct-Num'].startswith('too'):
                            graph.links.append(treatrow(row, source, target, smapping, tmapping))

        return graph        
        

    def add_link(self, sprot, spos, tprot, tpos, score, bidirectional = False):
        self.links.append({'Source-Start':spos[0],
                            'Source-End':spos[1],
                            'Source-Prot':sprot,
                            'Target-Start':spos[0],
                            'Target-End':spos[1],
                            'Target-Prot':sprot,
                            'Score':score
                            })
        if bidirectional:
            self.links.append({'Target-Start':spos[0],
                    'Target-End':spos[1],
                    'Target-Prot':sprot,
                    'Source-Start':spos[0],
                    'Source-End':spos[1],
                    'Source-Prot':sprot,
                    'Score':score
                    })


    def add_gene(self, name, length, color):
        self.genome.append({'name':name, 'color':color, 
                            'length':length})

    def genes_from_file(self, fname = 'Code/hiv_genome.txt'):
        with open(fname) as handle:
            for line in handle:
                if line.startswith('chr - '):
                    parts = line.strip().split()
                    self.add_gene(parts[2], int(parts[5]), parts[-1])               


    def make_figure(self, image_path, link_filter = None, link_key = None, 
                    link_max = None, SCRATCH = '/tmp/'):
        
        image_path = os.path.abspath(image_path)
        image_direc = image_path.rsplit(os.sep,1)[0]
        image_name = image_path.rsplit(os.sep,1)[1]
        print image_direc, image_name
        defaults = {'kayrotype':'hiv_genome.txt',
                    'outdir':image_direc,
                    'outfile':image_name,
                    'infile':'links.txt',
                    'radius':5000,
                    }

        circos_template = 'Code/circos.tem'
        copy_files = ('hiv_genome.txt', 'ideogram.conf',
                      'ticks.conf')
        tmp_loc = tempfile.mkdtemp(dir = SCRATCH)
        for f in copy_files:
            shutil.copy(os.path.join('Code', f),
                        os.path.join(tmp_loc, f))
        circos_file = os.path.join(tmp_loc, 'circos.conf')
        link_file = os.path.join(tmp_loc, 'links.txt')

        with open(link_file, 'w') as handle:
            link_iter = islice(ifilter(link_filter, sorted(self.links, key = link_key, reverse = True), 
                                ), link_max)
            
            for i, link in enumerate(link_iter):
                tdict = dict(**link)
                tdict['lnum'] = i
                handle.write('hivlink%(lnum)i %(Source-Prot)s %(Source-Start)i %(Source-End)i\n' % tdict)
                handle.write('hivlink%(lnum)i %(Target-Prot)s %(Target-Start)i %(Target-End)i\n' % tdict)

        with open(circos_template) as handle:
            ctemp = handle.read()
        
        rules = ''
        for ch in self.genome:
            rules += '<rule>\n'
            rules += 'condition = _CHR2_ eq "%s"\n' % ch['name']
            rules += 'color = %s\n' % ch['color']
            rules += '</rule>\n'
        defaults['rules'] = rules
    
        with open(circos_file, 'w') as handle:
            handle.write(ctemp % defaults)
        
        cmd = 'circos -conf circos.conf'
        with pushd(tmp_loc):
            args = shlex.split(cmd)
            call(args)
        shutil.rmtree(tmp_loc)


def guessing_figures(load_dir):
    
    min_rows = 25
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
        
    dist_cuts = numpy.arange(0, 60, 0.5)    
    link_cuts = numpy.arange(0.5, 1, 0.1)
    proteins.add(None)
    proteins = list(proteins) #so the order is consistent
    pnums = dict(zip(proteins, range(len(proteins))))
    vals = numpy.zeros((len(rows), 3))
    for i, row in enumerate(rows):
        vals[i-1,0] = pnums[row['Protein']]
        vals[i-1,1] = row['Linkage']
        vals[i-1,2] = row['3dDist']
    print vals

    for lval in link_cuts:
        

        fig, ax = pyplot.subplots(1)        
        for i, prot in enumerate(proteins):
            pvals = numpy.zeros_like(dist_cuts)
            if prot is None:
                tvals = vals
                prot = 'All'
            else:
                tvals = vals[vals[:,0] == (i-1), :]
            print tvals, tvals.size
            N_items = tvals.size
            for ind, dcut in enumerate(dist_cuts):
                X_successes = sum((tvals[:, 1] >= lval)*(tvals[:,2] <= dcut))
                n_draws = sum(tvals[:,2] <= dcut)                
                m_marked = sum(tvals[:, 1] >= lval)
                prob = hypergeo_cdf(X_successes, n_draws, m_marked, N_items)
            
                try:
                    pvals[ind-1] = -log(1-prob, 10)
                except:
                    print prob, ind, X_successes, n_draws, m_marked, N_items
                    #raise KeyError
                    pvals[ind-1] = numpy.nan
            ax.plot(dist_cuts[numpy.invert(numpy.isnan(pvals))], 
                    pvals[numpy.invert(numpy.isnan(pvals))], label = prot)
        pyplot.legend()
        fig.savefig(os.path.join(load_dir, 'figures', 'DvsLink-%i.png' % (10*lval,)))


