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

DATA_DIR = 'Data'
OUT_DIR = 'Results'
SEQ_LOC = os.path.join(DATA_DIR, 'ListFiles', 'sequences.list')
DOWNLOAD_XML = True

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
            part = f.split('.')[0]
            yield (os.path.join(load_dir, f), os.path.join(dump_dir, part + '.gi'))

@ruffus.files(seq_gen)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'RawSequences')),'get_sequence_xml')
def get_sequences(in_file, out_file):
    
    if not os.path.exists(out_file):
        with open(in_file) as handle:
            soup = BeautifulStoneSoup(handle.read())
        for seq, gi in extract_sequences(soup, XML = False, seq_only = True):
            with open(out_file, 'w') as handle:
                handle.write('>%s\n%s\n' % (gi, seq))


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


def seq_gen_fun():
    
    dump_dir = os.path.join(DATA_DIR, 'AASeqs') 
    mapping_file = os.path.join(DATA_DIR, 'ListFiles', 'mapping.txt')
    mapping_dict = make_mapping_dict(mapping_file)
    valid_ext = set(mapping_dict.values()) - set([None])
    for count, f in enumerate(xml_file_gen()):
        if count % 500 == 0:
            print 'Extracting:', count        
        base = f.split('.')[0]
        yield (f, mapping_file), [base + '.' + ext for ext in valid_ext], mapping_dict
            
@ruffus.files(seq_gen_fun)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'AASeqs')), 'get_sequence_xml')
def write_protein_sequences(in_files, out_files, mapping_dict):
    
    base = in_files[0].split('.')[0]
    for outdict in extract_features(in_files[0], mapping = mapping_dict):
        with open(base + '.' + outdict['name'], 'w') as handle:
            handle.write('>%s\n%s' % (base+'_'+outdict['name'], outdict['AA']))
    
@ruffus.files_re(os.path.join(DATA_DIR, 'KnownGenomes', '*'), 
                os.path.join(DATA_DIR, 'SubtypeBLAST', 'knownsubtypes.*'))
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'SubtypeBLAST')), 'write_protein_sequences')
def make_subtype_blast_db(in_files, out_files):

    known = {}
    with open(os.path.join(DATA_DIR, 'KnownGenomes', 'known.list')) as handle:
        for row in csv.DictReader(handle, delimiter = ','):
            known[row['gi']] = row['subtype']
    
    with open(os.path.join(DATA_DIR, 'SubtypeBLAST','knownsubtypes.fasta', 'w') as handle:
        for f in in_files:
            if f.endswith('.xml'):
                gi = f.rsplit(os.sep,1)[1].split('.')[0]
                for seq, _ in extract_sequences(soup, XML = False, seq_only = True):
                    handle.write('>%s_%s\n%s' % (gi, known[gi], seq))
    with pushd(os.path.join(DATA_DIR, 'SubtypeBLAST')):
        sh('formatdb -i knownsubtypes.fasta -p F')


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Code')
    parser.add_argument('--fresh', dest = 'fresh', action = 'store_true',
                         default = False)
    parser.add_argument('--no-seq', dest = 'noseqs', action = 'store_true',
                        default = False)
    parser.add_argument('--make-mapping', dest = 'makemapping', action = 'store_true',
                        default = False)
    parser.add_argument('--workers', dest = 'workers', default = 2,
                        action = 'store', type = int)
    args = parser.parse_args()
    
    if args.noseqs:
        DOWNLOAD_XML = False

    if args.fresh:
        touch_data()
        
    if args.makemapping:
        ruffus.pipeline_run([make_mappings])
    else:
        ruffus.pipeline_run([top_function], multiprocess = args.workers)



