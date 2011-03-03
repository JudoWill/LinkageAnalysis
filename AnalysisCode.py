import csv, os, os.path
import ruffus
from ruffus.proxy_logger import *
import urllib2, re, logging
from datetime import datetime
from itertools import *
from BeautifulSoup import BeautifulStoneSoup
import argparse
from Code.NCBIUtils import *
from Code.GeneralUtils import *
from Code.AlignUtils import *
from Code.ThreeDUtils import *
from Code.FigureTools import *
from functools import partial
from subprocess import call, check_call
import shlex, tempfile, shutil
from multiprocessing import Pool
from operator import itemgetter
from collections import defaultdict

DATA_DIR = 'Data'
OUT_DIR = 'Results'
FORCE_NEW = False
POOL_WORKERS = 3
WIN_SIZE = 50
WIN_OVERLAP = 25
LOG_BASE = os.path.join('logs', 'pipeline.log')
WIDTHS = range(1,5)
SUBTYPE = None
MIN_SEQS = 20
MIN_OVERLAP = 20
SPECIES_LIST = None
SPECIES_FILE = None


def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

@ruffus.files(SPECIES_FILE, SPECIES_FILE + '.sen')
def make_dirs(ifile, ofile):
    """Make all directories in the species list."""
    
    for species in SPECIES_LIST:
        for field, val in species.keys():
            if field.endswith('Dir'):
                safe_mkdir(val)
    touch(ofile)    

              
def align_gen():
    """Takes the generates the files which need to be aligned."""

    for species in SPECIES_LIST:
        if 'SequenceDir' in species and 'AlignmentDir' in species:
            seqdir = partial(os.path.join, species['SequenceDir'])
            aligndir = partial(os.path.join, species['AlignmentDir'])
            for f in os.listdir(species['SequenceDir']):
                if f.endswith('.fasta'):
                    name = f.split('.')[0]
                    ifile = seqdir(f)
                    ofiles = [aligndir(name+'.aln.fasta'), 
                                aligndir(name+'.aln')]
                    yield ifile, ofiles


@ruffus.files(align_gen)
@ruffus.follows(make_dirs)
def make_alignments(in_file, out_files, name):
    run_muscle(in_file, out_files[0])
    fasta2aln(out_files[0], out_files[1])

def lanl_align_pairs():
    load_dir = os.path.join(DATA_DIR, 'LANLSequences', 'Alignments')
    dump_dir = os.path.join(DATA_DIR, 'LinkageResults')
    aligns_present = [x.split('.')[0] for x in os.listdir(load_dir) if x.endswith('.aln')]
    ids_present = defaultdict(set)    
    for p1 in aligns_present:
        a1 = os.path.join(load_dir, p1+'.aln')
        aln = Alignment.alignment_from_file(a1)
        ids_present[p1] = set(aln.seqs.keys())
        if len(set(aln.seqs.values())) == 1:
            ids_present[p1] = set()
    aligns_present = [x for x in aligns_present if len(ids_present[x]) >= MIN_SEQS]
    

    for p1, p2 in product(sorted(aligns_present), repeat = 2):
    #for p1, p2 in zip(sorted(aligns_present), sorted(aligns_present)):
        if len(ids_present[p1] & ids_present[p2]) >= MIN_OVERLAP:

            a1 = os.path.join(load_dir, p1+'.aln')
            a2 = os.path.join(load_dir, p2+'.aln')
            
            d = os.path.join(dump_dir, p1+'--'+p2+'.res')
            s = os.path.join(dump_dir, p1+'--'+p2+'.sen')
            
            yield (a1, a2), (d, s)


@ruffus.files(lanl_align_pairs)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'LinkageResults')), 'make_alignments')
def calculate_lanl_linkages(in_files, out_files):
    print WIDTHS
    PredictionAnalysis(in_files[0], in_files[1], out_files[0], 
                        same = in_files[0] == in_files[1],
                        widths = WIDTHS)
    touch(out_files[1])

def scatter_files():
    struct_dir = os.path.join(DATA_DIR, 'ProteinStructures')
    linkage_dir = os.path.join(DATA_DIR, 'LinkageResults')
    align_dir = os.path.join(DATA_DIR, 'LANLSequences', 'Alignments')
    odir = os.path.join(DATA_DIR, 'ScatterResults')

    with open(os.path.join(struct_dir, 'mapping.txt')) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            struct_file = os.path.join(struct_dir, row['Structure'] + '.pdb')
            link_file = os.path.join(linkage_dir, row['Protein']+'--'+row['Protein'] + '.res')
            align_file = os.path.join(align_dir, row['Protein']+'.aln')
            if not os.path.exists(link_file) or not os.path.exists(align_file) \
                or not os.path.exists(struct_file):
                continue
            
            out_file = os.path.join(odir, row['Protein']+'--' + row['Structure'] + '.res')
            
            yield [struct_file, link_file, align_file], [out_file], row['Chain']



@ruffus.files(scatter_files)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'ScatterResults')),
                'calculate_lanl_linkages')
def generate_scatter(in_files, out_files, chain):
    args = in_files+out_files+[chain]
    create_scatter(*args)
    
    
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'ScatterResults', 'figures')),
                'generate_scatter')
def slice_scatters():
    
    guessing_figures(os.path.join(DATA_DIR, 'ScatterResults'))

def circos_files():
    link_path = os.path.join(DATA_DIR, 'LinkageResults')
    align_path = os.path.join(DATA_DIR, 'LANLSequences', 'Alignments')

    ifiles = []
    for f in os.listdir(link_path):
        if f.endswith('.res'):
            ifiles.append(os.path.join(link_path, f))
    for f in os.listdir(align_path):
        if f.endswith('.aln'):
            ifiles.append(os.path.join(align_path, f))
    yield ifiles, None


@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'CircosFigs')))
@ruffus.files(circos_files)
def circos_figs(ifile, ofile):
    dump_path = os.path.join(DATA_DIR, 'CircosFigs')
    load_path = os.path.join(DATA_DIR, 'LinkageResults')
    align_path = os.path.join(DATA_DIR, 'LANLSequences', 'Alignments')
    prots = set()
    for f in os.listdir(align_path):
        if f.endswith(f):
            prots.add(f)
    prots.add('All')

    lcuts = [0.5,0.6,0.7,0.8,0.9]
    graph = CircosGraph.load_from_dir(load_path, align_path)
    sort_fun = lambda x: x['Score']
    for lcut, prot in product(lcuts, prots):
        if prot == 'All':
            filter_fun = lambda x: x['Score'] >= lcut
        else:
            filter_fun = lambda x: x['Score'] >= lcut and x['Source-Prot'] == prot            
        path = os.path.join(dump_path, prot + '-%i.png' % (10*lcut,))
        print path
        graph.make_figure(path, link_key = sort_fun, link_filter = filter_fun)
    
    
    
        

  

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Code')
    parser.add_argument('--fresh', dest = 'fresh', action = 'store_true',
                         default = False)
    parser.add_argument('--make-mapping', dest = 'makemapping', action = 'store_true',
                        default = False)
    parser.add_argument('--workers', dest = 'workers', default = 3,
                        action = 'store', type = int)
    parser.add_argument('--max-width', dest = 'maxwidth', default = 1, action = 'store',
                        type = int)
    parser.add_argument('--data-dir', dest = 'datadir', default = 'OtherData', action = 'store')
    parser.add_argument('--processing-file', dest = 'processfile', default = None, action = 'store')
    parser.add_argument('--quiet', dest = 'quiet', action = 'store_true', default = False)
    parser.add_argument('--parse-align', dest = 'parsealign', action = 'store_true',
                        default = False)
    parser.add_argument('--link-lanl', dest = 'lanl', action = 'store_true', default = False)
    parser.add_argument('--filter-lanl', dest = 'filterlanl', action = 'store_true', default = False)
    parser.add_argument('--align', dest = 'alignments', action = 'store_true', default = False)
    parser.add_argument('--scatter-lanl', dest = 'lanlscatter', action = 'store_true', default = False)
    parser.add_argument('--circos-lanl', dest = 'lanlcircos', action = 'store_true', default = False)




    args = parser.parse_args()
    
    
    my_ruffus_logger = logging.getLogger('My_Ruffus_logger')
    my_ruffus_logger.setLevel(logging.INFO)
    
    fhandler = logging.handlers.RotatingFileHandler(LOG_BASE)
    my_ruffus_logger.addHandler(fhandler)
    
    if not args.quiet:
        shandler = logging.StreamHandler()
        my_ruffus_logger.addHandler(shandler)
    
    
    DATA_DIR = args.datadir

    if args.fresh:
        touch_data()
        FORCE_NEW = True
    WIDTHS = range(1, args.maxwidth+1)
    print WIDTHS

    if args.lanlscatter:
        ruffus.pipeline_run([slice_scatters], logger = my_ruffus_logger, multiprocess = args.workers)
    elif args.lanlcircos:
        ruffus.pipeline_run([circos_figs], logger = my_ruffus_logger)        
    elif args.alignments:
        ruffus.pipeline_run([make_alignments], logger = my_ruffus_logger, multiprocess = args.workers)        
    elif args.lanl:
        ruffus.pipeline_run([calculate_lanl_linkages], logger = my_ruffus_logger, multiprocess = args.workers)
    else:
        ruffus.pipeline_run([top_function], logger = my_ruffus_logger, multiprocess = args.workers)



