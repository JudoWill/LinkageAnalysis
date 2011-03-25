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

def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

@ruffus.files(SPECIES_FILE, SPECIES_FILE + '.sen')
def make_dirs(ifile, ofile):
    """Make all directories in the species list."""

    for species in SPECIES_LIST:
        for field, val in species.items():
            if field.endswith('Dir'):
                safe_mkdir(val)
    touch(ofile)

@ruffus.files(SPECIES_FILE, SPECIES_FILE + '.downloaded')
@ruffus.follows('make_dirs')
def download_data(ifile, ofile):    

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


@ruffus.jobs_limit(1)
@ruffus.files(align_gen)
@ruffus.follows('download_data', 'make_dirs')
def make_alignments(in_file, out_files):
    run_muscle(in_file, out_files[0])
    fasta2aln(out_files[0], out_files[1])

def align_pairs():
    """Yields pairs alignments for linkages."""

    def get_aligns(align_direc):
        """Gets the alignments present in the directory."""

        aligns = []
        for f in os.listdir(align_direc):
            if f.endswith('.aln'):
                aligns.append(f.split('.')[0])
        return aligns

    def get_ids_dict(align_direc):
        """Gets the IDS for each alignment in a directory."""
        
        aligndir = partial(os.path.join, align_direc)
        id_dict = defaultdict(set)
        for f in os.listdir(align_direc):
            if f.endswith('.aln'):
                aln = Alignment.alignment_from_file(aligndir(f))
                id_dict[f.split('.')[0]] = set(aln.seqs.keys())
        return id_dict
    


    for species in SPECIES_LIST:
        aligndir = partial(os.path.join, species['AlignmentDir'])
        linkagedir = partial(os.path.join, species['LinkageDir'])
        aligns_present = get_aligns(species['AlignmentDir'])
        align_ids = get_ids_dict(species['AlignmentDir'])
        widths = WIDTHS or species.get('WIDTHS', range(1,5))

        for p1, p2 in product(sorted(aligns_present), repeat = 2):
            print p1, p2, len(align_ids[p1] & align_ids[p2]), species.get('OVERLAP', 5)
            if len(align_ids[p1] & align_ids[p2]) > species.get('OVERLAP', 5):
                a1 = aligndir(p1 + '.aln')
                a2 = aligndir(p2 + '.aln')

                d = linkagedir(p1 + '--' + p2 + '.res')
                s = linkagedir(p1 + '--' + p2 + '.sen')
            
                yield (a1, a2), (d, s), widths


def tree_splitting():
    
    
    for species in SPECIES_LIST:
        if 'TreeDir' in species:
            afun = partial(os.path.join, species['AlignmentDir'])
            dfun = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            numcols = species.get('AlignmentCols', 100)
            ofiles = []
            for ind in xrange(numstraps):
                safe_mkdir(dfun('tree%i' % ind))
                ofiles.append(dfun('tree%i' % ind, 'infile'))
            ifiles = [afun(x) for x in os.listdir(afun('')) if x.endswith('.aln')]
            yield ifiles, ofiles, numcols

@ruffus.files(tree_splitting)
@ruffus.follows('make_alignments')
def tree_split(ifiles, ofiles, numcols):
    
    bigaln = Alignment.alignment_from_file(ifiles[0])
    for f in ifiles[1:]:
        bigaln.append_alignment(Alignment.alignment_from_file(f))

    bigaln.entropy_filter(numcols)
    for aln, f in izip(bigaln.bootstrap_columns(len(ofiles)), ofiles):
        aln.write_phylip(f)


def tree_run():

    for species in SPECIES_LIST:
        if 'TreeDir' in species:
            dfun = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            for ind in xrange(numstraps):
                ifile = dfun('tree%i' % ind, 'infile')
                ofile = dfun('tree%i' % ind, 'outtree')
                sfile = dfun('tree%i' % ind, 'process.sen')
                direc = dfun('tree%i' % ind)
                yield ifile, (ofile, sfile), direc

@ruffus.files(align_pairs)
@ruffus.follows('tree_split')
def process_trees(ifile, ofiles, direc):
    
    run_proml(direc)
    touch(ofiles[1])

def finished_trees():

    for species in SPECIES_LIST:
        if 'TreeDir' in species:
            dfun = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            ifiles = []
            mergedtree = dfun('intree')
            osfile = dfun('intree.sen')                
            for ind in xrange(numstraps):
                ifiles.append(dfun('tree%i' % ind, 'outtree'))
                ifiles.append(dfun('tree%i' % ind, 'process.sen'))
            yield ifiles, (mergedtree, osfile)

@ruffus.files(finished_trees)
@ruffus.follows('process_trees')
def tree_merge(ifiles, ofiles):
    
    with open(ofiles[0], 'w') as ohandle:
        treefiles = [x for x in ifiles if x.endswith('outtree')]
        for f in treefiles:
            with open(f) as handle:
                ohandle.write(handle.read())
    touch(ofiles[1])
    
    
    
                


@ruffus.files(align_pairs)
@ruffus.follows('make_alignments')
def calculate_linkages(in_files, out_files, widths):
    print in_files
    PredictionAnalysis(in_files[0], in_files[1], out_files[0], 
                        same = in_files[0] == in_files[1],
                        widths = widths)
    touch(out_files[1])

def linkage_merge():
    
    for species in SPECIES_LIST:
        circosdir = partial(os.path.join, species['CircosDir'])
        linkagedir = partial(os.path.join, species['LinkageDir'])
        
        ifiles = [linkagedir(x) for x in os.listdir(linkage_dir(''))]
        ofiles = [circosdir('FullAggregatedData.txt'),
                    circosdir('ShortAggregatedData.txt')]
        yield ifiles, ofiles

@ruffus.files(linkage_merge)
@ruffus.follows('calculate_linkages')
def merge_linkages(infiles, ofiles):
    
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
    parser.add_argument('--align', dest = 'alignments', action = 'store_true', default = False)
    parser.add_argument('--compare', dest = 'compare', action = 'store_true', default = False)
    parser.add_argument('--scatter-lanl', dest = 'lanlscatter', action = 'store_true', default = False)
    parser.add_argument('--circos-lanl', dest = 'lanlcircos', action = 'store_true', default = False)




    args = parser.parse_args()
    SPECIES_FILE = args.processfile

    with open(SPECIES_FILE) as handle:
        SPECIES_LIST = yaml.load(handle)


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
    elif args.link:
        ruffus.pipeline_run([merge_linkages], logger = my_ruffus_logger, multiprocess = args.workers)
    elif args.compare:
        ruffus.pipeline_run([compare_genomes], logger = my_ruffus_logger, multiprocess = args.workers)
    else:
        ruffus.pipeline_run([top_function], logger = my_ruffus_logger, multiprocess = args.workers)



