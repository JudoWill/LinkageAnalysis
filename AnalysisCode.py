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
from functools import partial
from subprocess import call, check_call
import shlex
from multiprocessing import Pool
from operator import itemgetter

DATA_DIR = 'Data'
OUT_DIR = 'Results'
SEQ_LOC = os.path.join('ListFiles', 'sequences.list')
DOWNLOAD_XML = True
FORCE_NEW = False
POOL_WORKERS = 3
WIN_SIZE = 50
WIN_OVERLAP = 25
LOG_BASE = os.path.join('logs', 'pipeline.log')

def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

def touch_data():
    for path, _, files in os.walk(DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            touch(os.path.join(path, f))

@ruffus.follows('calculate_linkages')
def top_function():
    pass

@ruffus.files(os.path.join('ListFiles', 'search_sentinal'),
              os.path.join('ListFiles', 'sequences.list'))
def get_sequence_ids(in_file, out_file):

    if not DOWNLOAD_XML:
        return

    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)
    logging.info('Retrieved %i sequences' % len(id_list))
    with open(out_file, 'w') as handle:
        handle.write('\n'.join(id_list))


@ruffus.files(os.path.join('ListFiles', 'sequences.list'),
              os.path.join(DATA_DIR, 'SequenceXML', 'download_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SequenceXML')), 
                'get_known_genotypes')
def get_sequence_xml(in_file, out_file):
    
    if not DOWNLOAD_XML:
        return
        
    dump_dir = os.path.join(DATA_DIR, 'SequenceXML')
    with open(in_file) as handle:
        seq_ids = set([x.strip() for x in handle])

    present = set([x.split('.')[0] for x in os.listdir(dump_dir)])
    dl_seqs = seq_ids - present
    for seq, gi in GetSeqs(dl_seqs, XML = True):
        logging.info('Retrieved %s from NCBI' % gi)
        with open(os.path.join(dump_dir, gi + '.xml'), 'w') as handle:
            handle.write('>%s\n%s\n' % (gi, seq))

    touch(out_file)

def seq_gen():
    load_dir = os.path.join(DATA_DIR, 'SequenceXML')
    dump_dir = os.path.join(DATA_DIR, 'RawSequences')
    for f in os.listdir(load_dir):
        if f.endswith('.xml'):
            yield os.path.join(load_dir, f)

@ruffus.files([os.path.join(DATA_DIR, 'KnownGenomes', 'download_sentinal'), 
                os.path.join(DATA_DIR, 'SequenceXML', 'download_sentinal')], 
                os.path.join(DATA_DIR, 'RawSequences', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'RawSequences')),'get_sequence_xml')
def get_sequences(in_files, out_file):

    dump_dir = os.path.join(DATA_DIR, 'RawSequences')    
    for f in seq_gen():
        gi = gi_from_path(f)
        out = os.path.join(dump_dir, gi + '.gi')
        if not os.path.exists(out) or FORCE_NEW:
            with open(f) as handle:
                soup = BeautifulStoneSoup(handle.read())      
            for seq, _ in extract_sequences(soup, XML = False, seq_only = True):
                logging.info('Extracting sequence from %s' % gi)
                with open(out, 'w') as handle:
                    handle.write('>%s\n%s\n' % (gi, seq.strip().upper()))
    touch(out_file)

@ruffus.files(os.path.join('ListFiles', 'known_subtypes.list'), 
                os.path.join(DATA_DIR, 'KnownGenomes', 'download_sentinal'))
@ruffus.follows(get_sequence_ids)
def get_known_genotypes(in_file, out):

    dump_dir = os.path.join('KnownGenomes')
    genomes = {}
    with open(in_file) as handle:
        for row in csv.DictReader(handle,
                                  fieldnames = ('ID', 'Subtype')):
            genomes[row['ID']] = row['Subtype']

    for xml, gi in GetSeqs(genomes.keys(), XML = True):
        logging.info('Retrieved %s from NCBI' % gi)
        with open(os.path.join(DATA_DIR, 'KnownGenomes', gi+'.xml'), 'w') as handle:
            handle.write(xml)
    
    touch(out)


def xml_file_gen():
    dirs = ('KnownGenomes', 'SequenceXML')
    files = []
    for d in dirs:
        di = os.path.join(DATA_DIR, d)
        for f in os.listdir(di):
            if f.endswith('.xml'):
                yield os.path.join(di, f)


@ruffus.merge([os.path.join(DATA_DIR, 'KnownGenomes', 'download_sentinal'), 
                os.path.join(DATA_DIR, 'SequenceXML', 'download_sentinal')], 
                os.path.join('ListFiles', 'mapping.txt'))
@ruffus.follows('get_known_genotypes', 'get_sequence_xml')
def make_mappings(in_files, out_file):
    
    mapping_dict = {}
    
    if os.path.exists(out_file):
        mapping_dict = make_mapping_dict(out_file)        

    count = 0
    for f in xml_file_gen():
        count += 1
        if count % 500 == 0:
            print 'Mapping: %i' % count
        for outdict in extract_features(f):
            if outdict['name'].lower() not in mapping_dict:
                mapping_dict[outdict['name'].lower()] = None
    
    with open(out_file, 'w') as handle:
        handle.write('%s\t%s\n' % ('key', 'name'))
        for key, val in sorted(mapping_dict.items(), key = lambda x: x[0]):
            handle.write('%s\t%s\n' % (key, val))
            
@ruffus.files([os.path.join(DATA_DIR, 'KnownGenomes', 'download_sentinal'), 
                os.path.join(DATA_DIR, 'SequenceXML', 'download_sentinal'),
                os.path.join('ListFiles', 'mapping.txt')],
                os.path.join(DATA_DIR, 'AASeqs', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'AASeqs')), 'get_sequence_xml')
def write_protein_sequences(in_files, out_file):
    
    mapping_dict = make_mapping_dict(in_files[-1])
    mapping_fun = partial(mapping_func, mapping_dict)
    dump_dir = os.path.join(DATA_DIR, 'AASeqs')
    done = set()
    for f in os.listdir(dump_dir):
        done.add(gi_from_path(f))

    for count, f in enumerate(xml_file_gen()):
        gi = gi_from_path(f)
        logging.info('Getting AAseqs from %s' % gi)
        if count % 500 == 0:
            logging.debug('Extracted %i' % count)
        any_in = False
        if gi not in done or FORCE_NEW:
            for outdict in extract_features(f, mapping = mapping_fun):
                any_in = True                
                loc = os.path.join(dump_dir, gi + '.' + outdict['name'])        
                with open(loc, 'w') as handle:
                    handle.write('>%s\n%s' % (gi+'_'+outdict['name'], outdict['AA'].strip().upper()))
            if not any_in:
                touch(os.path.join(dump_dir, gi))

    touch(out_file)
    
@ruffus.files(os.path.join('ListFiles', 'known_subtypes.list'), 
                os.path.join(DATA_DIR, 'SubtypeBLAST', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SubtypeBLAST')), 'write_protein_sequences')
def make_subtype_blast_db(in_file, out_file):
    BLAST_TYPE = guess_blast_computer_type()
    
    known = {}
    with open(in_file) as handle:
        for row in csv.DictReader(handle, delimiter = ','):
            known[row['gi']] = row['subtype']
    
    with open(os.path.join(DATA_DIR, 'SubtypeBLAST','knownsubtypes.fasta'), 'w') as handle:
        for f in os.listdir(os.path.join(DATA_DIR, 'KnownGenomes')):
            if f.endswith('.xml'):
                gi = gi_from_path(f)
                with open(os.path.join(DATA_DIR, 'KnownGenomes', f)) as inhandle:
                    soup = BeautifulStoneSoup(inhandle.read())
                for seq, _ in extract_sequences(soup, XML = False, seq_only = True):
                    handle.write('>%s_%s\n%s\n\n' % (gi, known[gi], seq.strip().upper()))
    with pushd(os.path.join(DATA_DIR, 'SubtypeBLAST')):
        logging.debug('Creating subtype BLAST database')
        cmd = make_blast_cmd('formatdb', None, 'knownsubtypes.fasta', None, 
                            blast_type = BLAST_TYPE, dbtype = 'nuc')       
        check_call(shlex.split(cmd))
    touch(out_file)

@ruffus.files(os.path.join('ListFiles', 'known_subtypes.list'), 
                os.path.join(DATA_DIR, 'TranslateBLAST', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'TranslateBLAST')), 'write_protein_sequences')
def make_translate_blast_db(in_file, out_file):
    BLAST_TYPE = guess_blast_computer_type()
    mapping_dict = make_mapping_dict(os.path.join('ListFiles', 'mapping.txt'))
    mapping = partial(mapping_func, mapping_dict)

    known = {}
    with open(in_file) as handle:
        for row in csv.DictReader(handle, delimiter = ','):
            known[row['gi']] = row['subtype']
    
    with open(os.path.join(DATA_DIR, 'TranslateBLAST','knowntranslations.fasta'), 'w') as handle:
        for f in os.listdir(os.path.join(DATA_DIR, 'KnownGenomes')):
            if f.endswith('.xml'):
                gi = gi_from_path(f)
                for outdict in extract_features(os.path.join(DATA_DIR, 'KnownGenomes', f), 
                                                mapping = mapping):
                    seq = outdict['AA']
                    prot = outdict['name']
                    handle.write('>%s_%s_%s\n%s\n\n' % (gi, known[gi], prot, seq.strip().upper()))
    
    
    with pushd(os.path.join(DATA_DIR, 'TranslateBLAST')):
        logging.debug('Creating Translating BLAST database')
        cmd = make_blast_cmd('formatdb', None, 'knowntranslations.fasta', None, 
                            blast_type = BLAST_TYPE, dbtype = 'prot')       
        check_call(shlex.split(cmd))
    touch(out_file)

            
@ruffus.files([os.path.join(DATA_DIR, 'SubtypeBLAST', 'processing_sentinal'),
                os.path.join(DATA_DIR, 'RawSequences', 'processing_sentinal')],
                os.path.join(DATA_DIR, 'SubtypeReports', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SubtypeReports')), 'make_subtype_blast_db', 'get_sequences')
def make_subtype_reports(in_files, out_file):
    BLAST_TYPE = guess_blast_computer_type()
        
    dump_dir = os.path.join(DATA_DIR, 'SubtypeReports')
    load_dir = os.path.join(DATA_DIR, 'RawSequences')
    db_path = os.path.join(DATA_DIR, 'SubtypeBLAST', 'knownsubtypes.fasta')
    done = set()
    for f in os.listdir(dump_dir):
        if f.endswith('.xml'):
            done.add(gi_from_path(f))
    
    for ind, f in enumerate(os.listdir(load_dir)):
        if f.endswith('.gi'):
            if ind % 500 == 0:
                logging.debug('Made %i Subtype Reports' % ind)
            gi = gi_from_path(f)
            if gi not in done or FORCE_NEW:
                with open(os.path.join(dump_dir, gi + '.fasta'), 'w') as handle:
                    for _, seq in fasta_iter(os.path.join(load_dir, f)):
                        count = 1
                        for block in OverlappingIterator(seq, WIN_SIZE, WIN_OVERLAP):
                            if sum([x.upper() == 'N' for x in block]) < len(block)*0.5:
                                handle.write('>%s\n%s\n' % (gi + '_' + str(count), ''.join(block)))
                            count += len(block)
                cmd = make_blast_cmd('blastn', db_path, os.path.join(dump_dir, gi + '.fasta'),
                                        os.path.join(dump_dir, gi + '.xml'))
                args = shlex.split(cmd)                
                logging.info('Calling BLAST for %s' % gi)
                retcode = call(args)
                if retcode != 0:
                    logging.warning('Got bad error code when calling BLASTn on %s' % gi)
                    os.remove(os.path.join(dump_dir, gi + '.fasta'))
                    os.remove(os.path.join(dump_dir, gi + '.xml'))
                
    touch(out_file)

@ruffus.files(os.path.join(DATA_DIR, 'SubtypeReports', 'processing_sentinal'),
                os.path.join('ListFiles', 'subtype_mapping.txt'))
@ruffus.follows('make_subtype_reports')
def process_subtype_reports(in_file, out_file):
    
    def get_subtypes(load_dir, gis, num_workers):
        if num_workers > 1:
            pool = Pool(processes = num_workers)
            fiter = [os.path.join(load_dir, x + '.xml') for x in gis]
            for f, subtype in izip(fiter, pool.imap(determine_subtype_element, fiter, chunksize = num_workers*5)):
                gi = gi_from_path(f)
                yield gi, subtype
        else:
            for gi in gis:
                yield gi, determine_subtype_element(os.path.join(load_dir, gi + '.xml'))

    
    load_dir = os.path.join(DATA_DIR, 'SubtypeReports')
    done = set()
    if FORCE_NEW or not os.path.exists(out_file):
        with open(out_file, 'w') as handle:
            handle.write('%s\t%s\n' % ('gi', 'Subtype'))
    
    with open(out_file) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            done.add(row['gi'])
    available = set()    
    for f in os.listdir(load_dir):
        if f.endswith('.xml'):
            available.add(gi_from_path(f))
    
    need = available - done
                
    with open(out_file, 'a') as handle:
        writer = csv.DictWriter(handle, ('gi', 'Subtype'), delimiter = '\t')
        for gi, subtype in get_subtypes(load_dir, need, POOL_WORKERS):
            need.discard(gi)
            writer.writerow({
                'gi':gi, 
                'Subtype': str(subtype)
                })
            if len(need) % 500 == 0:
                print 'Reports remaining: %i' % len(need)


@ruffus.files([os.path.join('ListFiles', 'subtype_mapping.txt'),
                os.path.join(DATA_DIR, 'AASeqs', 'processing_sentinal')],
                os.path.join(DATA_DIR, 'AlignmentDir', 'processing_sentinal'))
@ruffus.follows('make_subtype_reports', 'make_translate_blast_db',
                ruffus.mkdir(os.path.join(DATA_DIR, 'AlignmentDir')))                             
def make_alignments(in_files, out_file):
    
    
    load_dir = os.path.join(DATA_DIR, 'AASeqs')
    dump_base = os.path.join(DATA_DIR, 'AlignmentDir')
    prot_mapping = mapping_dict = make_mapping_dict(os.path.join('ListFiles', 'mapping.txt'))
    sub_mapping = defaultdict(set)
    
    with open(in_files[0]) as handle:
        c = 0
        for row in csv.DictReader(handle, delimiter = '\t'):
            sub_mapping[row['Subtype']].add(row['gi'])
            c+= 1

    uni_prots = set(prot_mapping.values())
    uni_prots.discard(None)
    
    for sub, gi_list in sub_mapping.items():
        print sub, gi_list
        logging.debug('Aligning Subtype %s' % sub)
        safe_mkdir(os.path.join(dump_base, sub))
        for prot in uni_prots:
            logging.debug('Aligning %s of subtype %s' % (prot, sub))
                
            dump_dir = os.path.join(dump_base, sub, prot)
            safe_mkdir(dump_dir)
            base_name = os.path.join(dump_dir, '%s-%s' % (sub, prot))
            if not os.path.exists(base_name + '.fasta'):
                logging.debug('Creating a new alignment')
                filenames = []
                for gi in gi_list:
                    if os.path.exists(os.path.join(load_dir, gi+'.'+prot)):
                        filenames.append(os.path.join(load_dir, gi+'.'+prot))
                if len(filenames) < 2:
                    continue
                print sub, prot, len(filenames)
                run_clustalw(filenames, 
                            base_name + '.fasta', 
                            base_name + '.dnd', 
                            base_name + '.align')
            else:
                logging.debug('Folding in new sequences')
                done_set = set([x for x,y in fasta_iter(base_name + '.fasta')])
                
                filenames = []
                for gi in gi_list - done_set:
                    #print gi
                    if os.path.exists(os.path.join(load_dir, gi+'.'+prot)):
                        filenames.append(os.path.join(load_dir, gi+'.'+prot))
                print len(done_set), len(gi_list - done_set), len(filenames)
                if len(filenames) < 2:
                    continue
                print 'merging', sub, prot, len(filenames)
                run_clustalw(filenames, 
                            base_name + 't.fasta', 
                            base_name + 't.dnd', 
                            base_name + 't.align')
                join_alignments(base_name + '.align', base_name + 't.align')
                join_fasta(filenames, base_name + '.fasta', strip = True, mode = 'a')
                os.remove(base_name + 't.fasta')
                os.remove(base_name + 't.dnd')
                os.remove(base_name + 't.align')
    
            
    
    
    touch(out_file)

def align_gen():
    for root, dirs, files in os.walk(os.path.join(DATA_DIR, 'AlignmentDir')):
        for f in files:
            if f.endswith('.align'):
                ipath = os.path.join(root, f)
                opath = os.path.join(root, f.split('.')[0] + '.aln')
                yield ipath, opath

@ruffus.files(align_gen)
@ruffus.follows('make_alignments')
def convert_alignments(in_file, out_file):
    convert_alignment(in_file, out_file)

def align_pairs():
    load_dir = os.path.join(DATA_DIR, 'AlignmentDir')
    dump_dir = os.path.join(DATA_DIR, 'LinkageResults')
    aligns_present = []
    for root, _, files in os.walk(load_dir):
        for f in files:
            if f.endswith('.aln'):
                parts = f.split('.')[0].split('-',1)
                aligns_present.append((parts[0], parts[1]))
    
    for subtype, prots in groupby(sorted(aligns_present), itemgetter(0)):
        needed = list([y for x, y in prots])
        for p1, p2 in product(needed, repeat = 2):
            print p1, p2, subtype            
            a1 = os.path.join(load_dir, subtype, p1, subtype+'-'+p1+'.aln')
            a2 = os.path.join(load_dir, subtype, p2, subtype+'-'+p2+'.aln')
            
            d = os.path.join(dump_dir, subtype+'-'+p1+'-'+p2+'.res')
            
            yield (a1, a2), d

@ruffus.files(align_pairs)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'LinkageResults')), 'convert_alignments')
def calculate_linkages(in_files, out_file):
    
    PredictionAnalysis(in_files[0], in_files[1], out_file, 
                        same = in_files[0] == in_files[1])
    




@ruffus.files(os.path.join(DATA_DIR, 'SubtypeReports', 'processing_sentinal'),
                os.path.join(DATA_DIR, 'SubtypeReports', 'simplfying_sentinal'))
def sanitize_xml(in_file, out_file):
    
    def file_gen(load_dir):
        for f in os.listdir(load_dir):
            if f.endswith('.xml'):
                yield os.path.join(load_dir, f)
                
    load_dir = os.path.join(DATA_DIR, 'SubtypeReports')
    
    if POOL_WORKERS == 1:
        for i, f in enumerate(file_gen(load_dir)):
            simplify_xml(f)
            if i % 500 == 0:
                print 'simplified: ', i
    else:
        pool = Pool(processes = POOL_WORKERS)
        for i, v in enumerate(pool.imap(simplify_xml, file_gen(load_dir),
                                chunksize = 100)):
            if i % 500 == 0:
                print 'simplified: ', i
                
    touch(out_file)

@ruffus.files([os.path.join(DATA_DIR, 'SubtypeReports', 'processing_sentinal'),
                os.path.join(DATA_DIR, 'RawSequences', 'processing_sentinal')],
                os.path.join(DATA_DIR, 'RawSequences', 'location_sentinal'))
def check_genome_locations(in_files, out_file):
        
    subtype_dir = os.path.join(DATA_DIR, 'SubtypeReports')
    genome_dir = os.path.join(DATA_DIR, 'RawSequences')
    
    count = 0
    for f in os.listdir(subtype_dir):
        if f.endswith('.xml'):        
            count += 1
            if count % 500 == 0:
                print 'Located:', count
            gi = gi_from_path(f)
            gpath = os.path.join(genome_dir, gi + '.gi')
            spath = os.path.join(subtype_dir, gi + '.xml')
            if os.path.exists(gpath):
                guess_location(spath, gpath, write_out = True)

    touch(out_file)

def direct_gen():
    
    def grouper(path):
        parts = path.split(os.sep)
        return parts[-1].split('-')[0]
    
    alignments = []
    for root, dirs, files in os.walk(os.path.join(DATA_DIR, 'AlignmentDir')):
        for f in files:
            if f.endswith('aln'):
                alignments.append(os.path.join(root, f))
                
    for subtype, aligns in groupby(sorted(alignments, key = grouper), key = grouper):
        yield list(aligns), os.path.join('ListFiles', 'LinkageReports', subtype+'.tdl')
        
@ruffus.files(direct_gen)
@ruffus.follows(ruffus.mkdir(os.path.join('ListFiles', 'LinkageReports')))
def make_overlap_reports(in_files, out_file):
    
    def prot_from_path(path):
        fname = path.split(os.sep)[-1].split('.')[0]
        parts = fname.split('-',1)
        return parts[1]
        
    prots = defaultdict(set)
    overlap = defaultdict(set)
    for f in in_files:
        prot = prot_from_path(f)
        with open(f) as handle:
            for line in handle:
                prots[prot].add(line.split('\t')[0])
    
    for p1, p2 in product(prots.keys(), repeat = 2):
        overlap[(p1, p2)] = prots[p1] & prots[p2]
        
    with open(out_file, 'w') as handle:
        fnames = ['Name'] + sorted(prots.keys())
        handle.write('\t'.join(fnames)+'\n')
        writer = csv.DictWriter(handle, fieldnames = fnames,
                                delimiter = '\t')
        rows = []
        for p1 in sorted(prots.keys()):
            rows.append({'Name':p1})
            for p2 in sorted(prots.keys()):
                rows[-1][p2] = len(overlap[(p1, p2)])
        writer.writerows(rows)
                
            
    
    
    

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Code')
    parser.add_argument('--fresh', dest = 'fresh', action = 'store_true',
                         default = False)
    parser.add_argument('--no-seq', dest = 'noseqs', action = 'store_true',
                        default = False)
    parser.add_argument('--make-mapping', dest = 'makemapping', action = 'store_true',
                        default = False)
    parser.add_argument('--workers', dest = 'workers', default = 3,
                        action = 'store', type = int)
    parser.add_argument('--simplify-xml', dest = 'simplifyxml', default = False,
                        action = 'store_true')
    parser.add_argument('--guess-locations', dest = 'guessloc', action = 'store_true',
                        default = False)
    parser.add_argument('--quiet', dest = 'quiet', action = 'store_true', default = False)
    parser.add_argument('--overlap-reports', dest = 'overlapreports', action = 'store_true',
                        default = False)
    parser.add_argument('--parse-align', dest = 'parsealign', action = 'store_true',
                        default = False)
    args = parser.parse_args()
    
    
    my_ruffus_logger = logging.getLogger('My_Ruffus_logger')
    my_ruffus_logger.setLevel(logging.INFO)
    
    fhandler = logging.handlers.RotatingFileHandler(LOG_BASE)
    my_ruffus_logger.addHandler(fhandler)
    
    if not args.quiet:
        shandler = logging.StreamHandler()
        my_ruffus_logger.addHandler(shandler)
    
    
    if args.noseqs:
        DOWNLOAD_XML = False

    if args.fresh:
        touch_data()
        FORCE_NEW = True
        
    if args.simplifyxml:
        ruffus.pipeline_run([sanitize_xml], logger = my_ruffus_logger)
    elif args.makemapping:
        ruffus.pipeline_run([make_mappings], logger = my_ruffus_logger)
    elif args.guessloc:
        ruffus.pipeline_run([check_genome_locations], logger = my_ruffus_logger)
    elif args.overlapreports:
        ruffus.pipeline_run([make_overlap_reports], logger = my_ruffus_logger)
    else:
        ruffus.pipeline_run([top_function], logger = my_ruffus_logger)



