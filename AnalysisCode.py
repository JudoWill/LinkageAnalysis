import csv, os, os.path
import ruffus
from ruffus.proxy_logger import *
import urllib2, re, logging
from datetime import datetime
from itertools import *
from BeautifulSoup import BeautifulStoneSoup
import argparse, yaml
from Code.NCBIUtils import *
from Code.GeneralUtils import *
from Code.AlignUtils import *
from Code.ThreeDUtils import *
from Code.FigureTools import *
from Code.CrossCompUtils import *
from Code.TreeUtils import *
from Code.RunUtils import *
from SequenceDownload import *
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
SPECIES_FILE = 'HIVData/HIVProcessing.yaml'
FileGen = partial(FileIter, SPECIES_FILE)
TOUCH_ONLY = False

def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

@ruffus.files(SPECIES_FILE, SPECIES_FILE + '.sen')
def make_dirs(ifile, ofile):
    """Make all directories in the species list."""

    with open(ifile) as handle:
        SPECIES_LIST = yaml.load(handle)

    for species in SPECIES_LIST:
        for field, val in species.items():
            if field.endswith('Dir'):
                safe_mkdir(val)
    touch(ofile)

@ruffus.files(SPECIES_FILE, SPECIES_FILE + '.downloaded')
@ruffus.follows('make_dirs')
def download_data(ifile, ofile):    

    with open(ifile) as handle:
        SPECIES_LIST = yaml.load(handle)

    urls = ('ftp://ftp.ncbi.nih.gov/genomes/Bacteria/',
            'ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/',)

    term_dict = {}    
    for species in SPECIES_LIST:
        if species.get('DOWNLOAD', False):
            term_dict[species['SpeciesName']] = species['DownloadDir']
    print term_dict
    if term_dict:
        check_NCBI(term_dict, urls)
        for species in SPECIES_LIST:
            if species.get('DOWNLOAD', False):
                unzip_dir(species['DownloadDir'])
                process_directory(species['DownloadDir'], out_direc = species['SequenceDir'])

    touch(ofile)
    

@ruffus.jobs_limit(1)
@ruffus.files(partial(FileGen, 'alignments'))
@ruffus.follows('download_data', 'make_dirs')
def make_alignments(in_file, out_files):

    if TOUCH_ONLY:
        touch_existing(out_files)
        return

    run_muscle(in_file, out_files[0])
    fasta2aln(out_files[0], out_files[1])


@ruffus.files(partial(FileGen, 'tree_splitting'))
@ruffus.follows('make_alignments')
def tree_split(ifiles, ofiles, numcols):

    if TOUCH_ONLY:
        touch_existing(ofiles)
        return
    
    bigaln = Alignment.alignment_from_file(ifiles[0])
    for f in ifiles[1:]:
        bigaln.append_alignment(Alignment.alignment_from_file(f))

    bigaln.entropy_filter(numcols)
    for aln, f in izip(bigaln.bootstrap_columns(len(ofiles)), ofiles):
        aln.write_phylip(f)


@ruffus.files(partial(FileGen, 'tree_run'))
@ruffus.follows('tree_split')
def process_trees(ifile, ofile, direc):

    if TOUCH_ONLY:
        touch_existing(ofile)
        return

    run_phylip(direc, 'proml')
    touch(ofile[1])

@ruffus.files(partial(FileGen, 'tree_merge'))
@ruffus.follows('process_trees')
def tree_merge(ifiles, ofile):

    if TOUCH_ONLY:
        touch_existing(ofile)
        return
    
    with open(ofile, 'w') as ohandle:
        treefiles = [x for x in ifiles if x.endswith('outtree')]
        for f in treefiles:
            with open(f) as handle:
                ohandle.write(handle.read())
    

@ruffus.files(partial(FileGen, 'tree_cons'))
@ruffus.follows('tree_merge')
def cons_tree(ifile, ofile, direc):

    if TOUCH_ONLY:
        touch_existing(ofile)
        return

    run_phylip(direc, 'consense')

@ruffus.files(partial(FileGen, 'merging_sequences'))
@ruffus.follows('cons_tree')
def merging_sequences(ifiles, ofiles):

    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    merge_sequences(ifiles[0], ofile[0], ifiles[1])


@ruffus.files(partial(FileGen, 'align_pairs'))
@ruffus.follows('make_alignments', 'merging_sequences')
def calculate_linkages(in_files, out_files, widths):

    if TOUCH_ONLY:
        touch_existing(out_files)
        return

    print in_files
    PredictionAnalysis(in_files[0], in_files[1], out_files[0], 
                        same = in_files[0] == in_files[1],
                        widths = widths)
    touch(out_files[1])

@ruffus.files(partial(FileGen, 'linkage_merge'))
@ruffus.follows('calculate_linkages')
def merge_linkages(infiles, ofiles):

    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    AggregateLinkageData(infiles, ofiles[0], ofiles[1])

def linkage_summarize():
    
    orgs = []
    fnames = []
    
    for species in SPECIES_LIST:
        circosdir = partial(os.path.join, species['CircosDir'])
        orgs.append(species['SpeciesName'])
        fnames.append(circosdir('ShortAggregatedData.txt'))
    ofiles = [os.path.dirname(SPECIES_FILE) + os.sep + 'ShortOverlap.txt',
            os.path.dirname(SPECIES_FILE) + os.sep + 'LongOverlap.txt']    
    yield fnames, ofiles, orgs


@ruffus.files(linkage_summarize)
@ruffus.follows('merge_linkages')
def compare_genomes(infiles, outfiles, orgnames):
    
    if TOUCH_ONLY:
        touch_existing(outfiles)
        return

    compare_linkages(infiles, orgnames, outfiles)
    
    

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
@ruffus.follows('calculate_linkages')
def generate_scatter(in_files, out_files, chain):
    args = in_files+out_files+[chain]
    create_scatter(*args)
    
    
@ruffus.follows('generate_scatter')
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
    parser.add_argument('--link', dest = 'link', action = 'store_true', default = False)
    parser.add_argument('--tree', dest = 'tree', action = 'store_true', default = False)
    parser.add_argument('--align', dest = 'alignments', action = 'store_true', default = False)
    parser.add_argument('--compare', dest = 'compare', action = 'store_true', default = False)
    parser.add_argument('--scatter-lanl', dest = 'lanlscatter', action = 'store_true', default = False)
    parser.add_argument('--circos-lanl', dest = 'lanlcircos', action = 'store_true', default = False)


    args = parser.parse_args()
    SPECIES_FILE = args.processfile
    
    DATA_DIR = args.datadir

    WIDTHS = range(1, args.maxwidth+1)
    print WIDTHS

    FileGen = partial(FileIter, SPECIES_FILE)
    TOUCH_ONLY = args.fresh


    if args.lanlscatter:
        ruffus.pipeline_run([slice_scatters], multiprocess = args.workers)
    elif args.lanlcircos:
        ruffus.pipeline_run([circos_figs])        
    elif args.alignments:
        ruffus.pipeline_run([make_alignments], multiprocess = args.workers)        
    elif args.link:
        ruffus.pipeline_run([merge_linkages], multiprocess = args.workers)
    elif args.compare:
        ruffus.pipeline_run([compare_genomes], multiprocess = args.workers)
    elif args.tree:
        ruffus.pipeline_run([cons_tree], multiprocess = args.workers)
    else:
        ruffus.pipeline_run([top_function], multiprocess = args.workers)



