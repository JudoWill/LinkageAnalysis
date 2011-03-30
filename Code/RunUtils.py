import sqlite3
from types import ListType, TupleType, StringType
from itertools import product, chain, repeat
from functools import partial
from collections import defaultdict
from operator import contains, itemgetter
import os.path, os
import hashlib, yaml
from GeneralUtils import safe_mkdir, take
from AlignUtils import Alignment
from memorised.decorators import memorise

def hexhash(path):
    """Returns the md5 hexhash of the file at the path."""

    m = hashlib.md5()
    with open(path) as handle:
        for line in handle:
            m.update(line)
    return m.hexdigest()

def try_connect(fname, num_tries = 500, wait_length = 10):
    
    for x in xrange(num_tries):
        try:
            con = sqlite3.connect(fname, wait_length)
        except sqlite3.OperationalError:
            continue
        return con
    return sqlite3.connect(fname, wait_length)
            


def need_to_do(fname, ifiles, ofiles, *args, **kwargs):
    """Checks whether the files need to be updated."""

    ohashes = {}
    if type(ofiles) == ListType or type(ofiles) == TupleType:
        for f in ofiles:
            if not os.path.exists(f):
                return True, 'Missing file: %s' % f
            else:
                ohashes[f] = hexhash(f)
    else:
        if not os.path.exists(ofiles):
            return True, 'Missing file: %s' % f
        else:
            ohashes[ofiles] = hexhash(ofiles)

    if type(ifiles) == ListType or type(ifiles) == TupleType:
        ihashes = dict([(x, hexhash(x)) for x in ifiles])
    else:
        ihashes = dict()
        ihashes[ifiles] = hexhash(ifiles)

    con = try_connect('filedata.sql')
    stri = "select * from dep where fname=? and spath=? and shash=? and dpath=? and dhash=?"
    iterable = product(ihashes.iteritems(), ohashes.iteritems())
    for (spath, shash), (dpath, dhash) in iterable:
        rows = con.execute(stri, (fname, spath, shash, dpath, dhash))
        r = rows.fetchone()
        #print r
        if r is None:
            return True, 'Files out of date!'
    
    return False, 'All files up to date!'

def add_complete_files(fname, ifiles, ofiles):
    """A function for adding files to the database"""

    def del_gen(fname, ihashes, ohashes):
        for ip, op in product(ihashes.iterkeys(), ohashes.iterkeys()):
            yield (fname, ip, op)

    def in_gen(fname, ihashes, ohashes):
        for (ip, ih), (op, oh) in product(ihashes.iteritems(), ohashes.iteritems()):
            yield (fname, ip, ih, op, oh)

    ihashes = dict([(x, hexhash(x)) for x in ifiles])
    ohashes = dict([(x, hexhash(x)) for x in ofiles])

    con = try_connect('filedata.sql')
    dstr = "delete from dep where fname=? and spath=? and dpath=?"
    istr = "insert into dep values (?,?,?,?,?)"
    with con:
        con.executemany(dstr, del_gen(fname, ihashes, ohashes))
        con.executemany(istr, in_gen(fname, ihashes, ohashes))

    


def create_filedatabase():
    
    if not os.path.exists('filedata.sql'):
        con = sqlite3.connect('filedata.sql')
        con.execute('create table dep (fname text, spath text, shash text, dpath text, dhash text)')


def get_ids_dict(align_direc):
    """Gets the IDS for each alignment in a directory."""
    
    aligndir = partial(os.path.join, align_direc)
    id_dict = defaultdict(set)
    for f in os.listdir(align_direc):
        if f.endswith('.aln'):
            aln = Alignment.alignment_from_file(aligndir(f))
            id_dict[f.split('.')[0]] = set(aln.seqs.keys())
    return id_dict


def check_fields(species, fields):
    """Checks to make sure the fields are present"""

    for f in fields:
        if f not in species:
            return False
    return True

def FileIter(species_file, funcname):
    """A large iterator for all functions in AnalysisCode"""

    with open(species_file) as handle:
        SPECIES_LIST = yaml.load(handle)

    field_dict = {'alignments':('AlignmentDir', 'LinkageDir'),
                'align_pairs':('AlignmentDir', 'LinkageDir'),
                'tree_splitting':('TreeDir', 'AlignmentDir'),
                'tree_run':('TreeDir',),
                'tree_merge':('TreeDir',),
                'tree_cons':('TreeDir',),
                'linkage_merge':('LinkageDir', 'CircosDir'),
                'compare_genomes':('CircosDir', 'ComparingGenomes'),
                'merging_sequences': ('AlignmentDir', 'TreeDir', 'MergedDir'),
                'make_dirs':tuple(),
                'download_data':tuple(),
                'convert_linkages':('RefGenome', 'LinkageDir', 'AlignmentDir')
            }

    for species in SPECIES_LIST:
        if not check_fields(species, field_dict[funcname]):
            continue

        if funcname == 'make_dirs':
            yield species_file, species_file + '.sen'
            break #this only needs to be run once for each file!

        elif funcname == 'download_data':
            yield species_file, species_file + '.downloaded'
            break #this only needs to be run once for each file!

        elif funcname == 'alignments':
            seqdir = partial(os.path.join, species['SequenceDir'])
            aligndir = partial(os.path.join, species['AlignmentDir'])
            for f in sorted(os.listdir(seqdir(''))):
                if not f.endswith('skip'):
                    name = f.split('.')[0]
                    ifile = seqdir(f)
                    ofiles = [aligndir(name+'.aln.fasta'),
                                aligndir(name+'.aln')]
                    yield ifile, ofiles
        
        elif funcname == 'merging_sequences':
            aligndir = partial(os.path.join, species['AlignmentDir'])
            mergedir = partial(os.path.join, species['MergedDir'])
            refgenome = species.get('RefGenome', None)
            tfile = os.path.join(species['TreeDir'], 'outtree')
            files = sorted([x for x in os.listdir(aligndir('')) if x.endswith('.aln')])
            for f in files:
                yield [aligndir(f), tfile], [mergedir(f), mergedir(f+'.fasta')], [refgenome]

        elif funcname == 'align_pairs':
            if 'MergedDir' in species:
                aligndir = partial(os.path.join, species['MergedDir'])
            else:
                aligndir = partial(os.path.join, species['AlignmentDir'])
            linkdir = partial(os.path.join, species['LinkageDir'])

            aligns = [x.split('.')[0] for x in os.listdir(aligndir('')) if x.endswith('.aln')]
            align_ids = get_ids_dict(aligndir(''))
            overlap = species.get('OVERLAP', 5)
            widths = species.get('WIDTHS', range(1,5))
            for p1, p2 in product(sorted(aligns), repeat = 2):
                if len(align_ids[p1] & align_ids[p2]) >= overlap:
                    a1 = aligndir(p1 + '.aln')
                    a2 = aligndir(p2 + '.aln')
                    d = linkdir(p1 + '--' + p2 + '.res')
                    s = linkdir(p1 + '--' + p2 + '.sen')                

                    yield (a1, a2), (d, s), widths

        elif funcname == 'convert_linkages':

            ref_genome = species['RefGenome']
            for (a1, a2), (link_file, _), _ in FileIter(species_file, 'align_pairs'):
                out_file = link_file+'.conv'
                sen_file = link_file+'.conv.sen'
                if os.path.exists(link_file):
                    yield (a1, a2, link_file), (out_file, sen_file), ref_genome




        elif funcname == 'tree_splitting':
            aligndir = partial(os.path.join, species['AlignmentDir'])
            treedir = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            numcols = species.get('AlignmentCols', 100)
            ofiles = []

            for ind in xrange(numstraps):
                safe_mkdir(treedir('tree%i' % ind))
                ofiles.append(treedir('tree%i' % ind, 'infile'))
            ifiles = [aligndir(x) for x in os.listdir(aligndir('')) if x.endswith('.aln')]
            
            yield ifiles, ofiles, numcols

        elif funcname == 'tree_run':
            treedir = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            for ind in xrange(numstraps):
                ifile = treedir('tree%i' % ind, 'infile')
                ofile = treedir('tree%i' % ind, 'outtree')
                osfile = treedir('tree%i' % ind, 'outtree.sen')
                direc = treedir('tree%i' % ind)
                yield ifile, (ofile, osfile), direc

        elif funcname == 'tree_merge':
            treedir = partial(os.path.join, species['TreeDir'])
            numstraps = species.get('Bootstraps', 100)
            ifiles = [treedir('tree%i' % x, 'outtree') for x in xrange(numstraps)]
            mergedtree = treedir('intree')
            yield ifiles, mergedtree

        elif funcname == 'tree_cons':
            dfun = partial(os.path.join, species['TreeDir'])
            mergedtree = dfun('intree')
            otree = dfun('outtree')
            direc = dfun('')
            yield mergedtree, otree, direc
            

        elif funcname == 'linkage_merge':
            linkdir = partial(os.path.join, species['LinkageDir'])
            cirdir = partial(os.path.join, species['CircosDir'])
            
            ifiles = [linkdir(x) for x in os.listdir(linkdir('')) if x.endswith('.conv')]
            ofiles = [cirdir('FullAggregatedData.txt'),
                        cirdir('ShortAggregatedData.txt')]

            yield ifiles, ofiles



def mark_present_as_correct(species_file, need_delete = False, block_size = 10000):
    """A utility function for marking all present files as "correct" """

    def group_iter(species_file):
        fnames = ('alignments', 'align_pairs', 'tree_splitting', 'tree_run',
                   'tree_merge', 'linkage_merge')

        for name in fnames:
            for tup in FileIter(species_file, name):
                ifiles = tup[0]
                ofiles = tup[1]
                if type(ifiles) == StringType:
                    ifiles = [ifiles]
                if type(ofiles) == StringType:
                    ofiles = [ofiles]
                pfiles = [os.path.exists(x) for x in chain(ifiles, ofiles)]
                if all(pfiles):
                    print 'marking', ifiles, ofiles, name
                    for ifile, ofile in product(ifiles, ofiles):
                        yield name, ifile, hexhash(ifile), ofile, hexhash(ofile)

    con = sqlite3.connect('filedata.sql')
    dstr = "delete from dep where fname=? and spath=? and dpath=?"
    istr = "insert into dep values (?,?,?,?,?)"
    dgetter = itemgetter(0, 1, 3)
    iterable = group_iter(species_file)
    block = take(block_size, iterable)
    while block:
        
        if need_delete:
            con.executemany(dstr, map(dgetter, block))            

        con.executemany(istr, block)
        con.commit()
       
        block = take(block_size, iterable)






