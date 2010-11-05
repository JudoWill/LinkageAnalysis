from paver.easy import sh
import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime
from itertools import imap, islice
from suds.client import Client
from BeautifulSoup import BeautifulStoneSoup
import argparse
from Code.NCBIUtils import *

DATA_DIR = 'Data'
OUT_DIR = 'Results'


def touch_data():
    for path, _, files in os.walk(DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            sh('touch %s' % os.path.join(path, f))



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
@ruffus.follows(get_sequence_ids)
def get_sequences(in_file, out_file):

    dump_dir = os.path.join(DATA_DIR, 'RawSequences')
    with open(in_file) as handle:
        seq_ids = set([x.strip() for x in handle])

    present = set([x.split('.')[0] for x in os.listdir(dump_dir)])
    dl_seqs = seq_ids - present
    seq_iter = iter(dl_seqs)

    client = Client('http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/efetch_seq.wsdl',
                    retxml = True)
    block_size = 100
    block = take(block_size, seq_iter)
    count = 0
    while block:
        for seq, gi in GetSeqs(client, block):
            with open(os.path.join(dump_dir, gi + '.gi'), 'w') as handle:
                handle.write('>%s\n%s\n' % (gi, seq))
        count += len(block)
        print count, len(dl_seqs)
        block = take(block_size, seq_iter)

    sh('touch %s' % out_file)

######Utility Functions!







def take(N, iterable):
    return list(islice(iterable, N))


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Linkage Analysis Code')
    parser.add_argument('--fresh', dest = 'fresh', action = 'store_true',
                         default = False)
    args = parser.parse_args()


    if args.fresh:
        touch_data()

    ruffus.pipeline_run([top_function])


