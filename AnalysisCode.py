import csv, os, os.path, sys
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

    parser.add_argument('--do-all', dest = 'doall', action = 'store_true', default = False)


    args = parser.parse_args()
    SPECIES_FILE = args.processfile
    
    DATA_DIR = args.datadir

    WIDTHS = range(1, args.maxwidth+1)
    print WIDTHS

    FileGen = partial(FileIter, SPECIES_FILE)
    TOUCH_ONLY = args.fresh

    if args.doall:
        current_species_files = ('HIVData/HIVProcessing.yaml',
                                'HCVSeqs/HCVProcessing.yaml',
                                'BacterialData/BacterialProcessing.yaml',
                                'HIVData/HIVProcessing_long.yaml')
        
        for f in current_species_files:
            inargs = shlex.split('python '+' '.join(sys.argv) + ' --processing-file ' + f)
            inargs = [x for x in inargs if x != '--do-all']
            check_call(inargs)

        raise SystemExit

    



def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

@ruffus.files(partial(FileGen, 'make_dirs'))
def make_dirs(ifile, ofile):
    """Make all directories in the species list.
    
    Arguments:
    ifile -- Input yaml species file
    ofile -- Directory sentinal
    """

    with open(ifile) as handle:
        SPECIES_LIST = yaml.load(handle)

    for species in SPECIES_LIST:
        for field, val in species.items():
            if field.endswith('Dir'):
                safe_mkdir(val)
    touch(ofile)

@ruffus.files(partial(FileGen, 'download_data'))
@ruffus.follows('make_dirs')
def download_data(ifile, ofile):    
    """Downloads data from NCBI.
    
    Arguments:
    ifile -- Input yaml species file
    ofile -- Download sentinal
    """
    
    
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
def make_alignments(ifile, ofiles):
    """Align sequences.
    
    Uses the muscle program to create multiple alignments from input fasta 
    files.
    
    Arguments:
    ifile -- An input fasta file.
    ofiles -- A tuple of output files: (fasta-file, tab-delimited-file)
    """

    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    run_muscle(ifile, ofiles[0])
    if os.path.exists(ofiles[0]):
        fasta2aln(ofiles[0], ofiles[1])
    else:
        shutil.move(ifile, ifile+'.skip')


@ruffus.files(partial(FileGen, 'tree_splitting'))
@ruffus.follows('make_alignments')
def tree_split(ifiles, ofiles, numcols):
    """Joins alignments and then finds columns for tree generation.
    
    This implements the bootstraping algorithm for generating trees. All 
    input alignments are joined based on thier key-values and then the highest
    entropy columns are put into a single alignment. These are then selected 
    WITH REPLACEMENT and split into multiple files for parallel processing.
    
    Arguments:
    ifiles -- A list of tab-delimited input files to join.
    ofiles -- A list of input files to create ... this should be the number of 
                bootstraps to perform.
    numcols -- An integer indicating the number of columns to including in each
                bootstrap.
    """

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
@ruffus.jobs_limit(3)
def process_trees(ifile, ofile, direc):
    """Run phylip proml on a tree.
    
    Arguments:
    ifile -- The phylip-formatted alignment to process
    ofiles -- A tuple of files to create (tree-file, sentinal-file)
    direc -- The directory which contains the input files.
    """

    if TOUCH_ONLY:
        touch_existing(ofile)
        return

    run_phylip(direc, 'proml')
    touch(ofile[1])

@ruffus.files(partial(FileGen, 'tree_merge'))
@ruffus.follows('process_trees')
def tree_merge(ifiles, ofile):
    """Merge a set of trees into one file.
    
    This step is needed because the consensus program requires that all trees
    be in one file but in the previous step we split them to achieve parallel 
    execution.
    
    Arguments:
    ifiles -- A list of tree files to join.
    ofile -- The output file to create.
    """

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
    """Run the phylip consensus program.
    
    Arguments:
    ifile -- The input JOINED tree file.
    ofile -- The consensus tree produced.
    """
    

    if TOUCH_ONLY:
        touch_existing(ofile)
        return

    run_phylip(direc, 'consense')

@ruffus.files(partial(FileGen, 'merging_sequences'))
@ruffus.follows('cons_tree')
def merging_sequences(ifiles, ofiles, excluded):
    """Merges similar sequences based on an input tree.
    
    This is needed for the "phylogentic correction" step. This algorithm
    joins sequences which are determined to be similar by thier bootstrap 
    values in a phylogenetic tree. It creates a new alignment file.
    
    Arguments:
    ifiles -- A tuple (alignment-file, tree-file)
    ofiles -- A tuple (tab-delimited alingment-file, fasta-formated alignment)
    excluded -- A list of sequences to NEVER merge. This is usually just the 
                reference sequence since it will be needed in later steps.
    """

    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    merge_sequences(ifiles[0], ofiles[0], ifiles[1], excluded = excluded)


@ruffus.files(partial(FileGen, 'align_pairs'))
@ruffus.follows('make_alignments', 'merging_sequences')
def calculate_linkages(ifiles, ofiles, widths):
    """Calculate linkages based on a pair of input files.
    
    Arguments:
    ifiles -- A tuple of alignment files to compare.
    ofiles -- A tuple (output-linkages, sentinal-file)
    widths -- A list of widths to check ... traditionally this is [1,2,3,4,5]
    """


    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    print ifiles,  min(widths, WIDTHS, key = len)
    PredictionAnalysis(ifiles[0], ifiles[1], ofiles[0], 
                        same = ifiles[0] == ifiles[1],
                        widths = min(widths, WIDTHS, key = len))
    touch(ofiles[1])



@ruffus.files(partial(FileGen, 'convert_linkages'))
@ruffus.follows('calculate_linkages')
def fix_numbering(ifiles, ofiles, ref_genome):
    """Converts the numbering from "alignment-space" to "sequence-space"
    
    Since the numbering in the linkage results are dependant on the exact 
    MSA used it needs to be converted into something that is more consistent.
    
    Arguments:
    ifiles -- A tuple of input files:
            (Source Alignment, Target Aligment, linkage-file)
    ofiles -- A tuple of output files:
            (Output-linkage-file, sentinal-file)
    """
    if TOUCH_ONLY:
        touch_existing(ofiles)
        return
    
    convert_numbering(ifiles[0], ifiles[1], ifiles[2], ofiles[0], ref_genome)
    touch(ofiles[1])


@ruffus.files(partial(FileGen, 'linkage_merge'))
@ruffus.follows('fix_numbering')
def merge_linkages(infiles, ofiles):
    """Merges linkage files into one summary file.
    
    Arguments:
    infiles -- A list of linkage files to merge.
    ofiles -- A tuple of output files:
        (Full-Linkage, Short-Linkage)
    """

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
def compare_genomes(ifiles, ofiles, orgnames):
    
    if TOUCH_ONLY:
        touch_existing(ofiles)
        return

    compare_linkages(ifiles, orgnames, ofiles)
    
    

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




    
        

  




