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
              os.path.join(DATA_DIR, 'RawSequences', 'download_sentinal'))
@ruffus.follows('get_known_genotypes')
def get_sequences(in_file, out_file):

    dump_dir = os.path.join(DATA_DIR, 'RawSequences')
    with open(in_file) as handle:
        seq_ids = set([x.strip() for x in handle])

    present = set([x.split('.')[0] for x in os.listdir(dump_dir)])
    dl_seqs = seq_ids - present
    for seq, gi in GetSeqs(dl_seqs):
        with open(os.path.join(dump_dir, gi + '.gi'), 'w') as handle:
            handle.write('>%s\n%s\n' % (gi, seq))

    sh('touch %s' % out_file)

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
    args = parser.parse_args()


    if args.fresh:
        touch_data()

    ruffus.pipeline_run([top_function])



