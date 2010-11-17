import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime
from itertools import imap, islice, izip
from BeautifulSoup import BeautifulStoneSoup
import argparse
from Code.NCBIUtils import *
from Code.GeneralUtils import *
from Code.AlignUtils import *
from functools import partial
from subprocess import call, check_call
import shlex
from multiprocessing import Pool


DATA_DIR = 'Data'
OUT_DIR = 'Results'
SEQ_LOC = os.path.join('ListFiles', 'sequences.list')
DOWNLOAD_XML = True
FORCE_NEW = False
POOL_WORKERS = 3
WIN_SIZE = 50
WIN_OVERLAP = 25


def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

def touch_data():
    for path, _, files in os.walk(DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            touch(os.path.join(path, f))

@ruffus.follows('make_alignments')
def top_function():
    pass

@ruffus.files(os.path.join('ListFiles', 'search_sentinal'),
              os.path.join('ListFiles', 'sequences.list'))
def get_sequence_ids(in_file, out_file):



    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)
    print len(id_list)
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
                with open(out, 'w') as handle:
                    handle.write('>%s\n%s\n' % (gi, seq.strip().upper()))
    touch(out_file)

@ruffus.files(os.path.join(DATA_DIR, 'KnownGenomes', 'known.list'), 
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
        if count % 500 == 0:
            print 'Extracting:', count
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
    
@ruffus.files(os.path.join(DATA_DIR, 'KnownGenomes', 'known.list'), 
                os.path.join(DATA_DIR, 'SubtypeBLAST', 'processing_sentinal'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SubtypeBLAST')), 'write_protein_sequences')
def make_subtype_blast_db(in_file, out_file):
    BLAST_TYPE = guess_blast_computer_type()
    
    known = {}
    with open(os.path.join(DATA_DIR, 'KnownGenomes', 'known.list')) as handle:
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
        cmd = make_blast_cmd('formatdb', None, 'knownsubtypes.fasta', None, 
                            blast_type = BLAST_TYPE, dbtype = 'nuc')       
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
                print ind
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
                                
                retcode = call(args)
                if retcode != 0:
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
@ruffus.follows('make_subtype_reports', 
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
        #print sub, gi_list
        safe_mkdir(os.path.join(dump_base, sub))
        for prot in uni_prots:
            
                
            dump_dir = os.path.join(dump_base, sub, prot)
            safe_mkdir(dump_dir)
            base_name = os.path.join(dump_dir, '%s-%s' % (sub, prot))
            if not os.path.exists(base_name + '.fasta'):
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
    args = parser.parse_args()
    
    if args.noseqs:
        DOWNLOAD_XML = False

    if args.fresh:
        touch_data()
        FORCE_NEW = True
        
    if args.makemapping:
        ruffus.pipeline_run([make_mappings])
    else:
        ruffus.pipeline_run([top_function])



