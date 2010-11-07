from paver.easy import sh
import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime
from itertools import imap, islice
from BeautifulSoup import BeautifulStoneSoup
import argparse
from Code.NCBIUtils import *

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

@ruffus.follows('get_sequences')
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
            yield (os.path.join(load_dir, f), os.path.join(dump_dir, part + 'gi'))

@ruffus.files(seq_gen)
@ruffus.follows(ruffus.mkdir(os.path.join(DATA_DIR, 'RawSequences')),'get_sequence_xml')
def get_sequences(in_file, out_file):
    
    with open(in_file) as handle:
        soup = BeautifulStoneSoup(handle.read())
    for seq, gi in extract_sequences(soup):
        with open(out_file, 'w') as handle:
            handle.write('>%s\n%s\n' % (gi, seq))


@ruffus.files(os.path.join(DATA_DIR, 'KnownGenomes', 'known.list'), None)
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








if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Code')
    parser.add_argument('--fresh', dest = 'fresh', action = 'store_true',
                         default = False)
    parser.add_argument('--no-seq', dest = 'noseqs', action = 'store_true',
                        default = False)
    args = parser.parse_args()
    
    if args.noseqs:
        DOWNLOAD_XML = False

    if args.fresh:
        touch_data()

    ruffus.pipeline_run([top_function])



