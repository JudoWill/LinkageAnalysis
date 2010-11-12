from paver.easy import sh, pushd
import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime
from itertools import imap, islice
from BeautifulSoup import BeautifulStoneSoup
import argparse
from Code.NCBIUtils import *
from Code.GeneralUtils import *
from functools import partial

DATA_DIR = 'Data'
OUT_DIR = 'Results'
SEQ_LOC = os.path.join(DATA_DIR, 'ListFiles', 'sequences.list')
DOWNLOAD_XML = True
FORCE_NEW = False
POOL_WORKERS = 1

def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

def touch_data():
    for path, _, files in os.walk(DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            touch(os.path.join(path, f))

@ruffus.follows('make_subtype_blast_db')
def top_function():
    pass

@ruffus.files(os.path.join(DATA_DIR, 'ListFiles', 'search_sentinal'),
              os.path.join(DATA_DIR, 'ListFiles', 'sequences.list'))
def get_sequence_ids(in_file, out_file):



    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)
    print len(id_list)
    with open(out_file, 'w') as handle:
        handle.write('\n'.join(id_list))


@ruffus.files(os.path.join(DATA_DIR, 'ListFiles', 'sequences.list'),
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
                os.path.join(DATA_DIR, 'ListFiles', 'mapping.txt'))
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
                os.path.join(DATA_DIR, 'ListFiles', 'mapping.txt')],
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
                os.path.join(DATA_DIR, 'SubtypeBLAST', 'knownsubtypes.*'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SubtypeBLAST')), 'write_protein_sequences')
def make_subtype_blast_db(in_file, out_files):

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
                    handle.write('>%s_%s\n%s' % (gi, known[gi], seq.strip().upper()))
    with pushd(os.path.join(DATA_DIR, 'SubtypeBLAST')):
        try:
            sh('formatdb -i knownsubtypes.fasta -p F')
        except:
            sh('makeblastdb -in knownsubtypes.fasta -dbtype nucl')


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



