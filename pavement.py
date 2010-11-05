from paver.easy import *
import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime
from itertools import imap, islice
from suds.client import Client
from BeautifulSoup import BeautifulStoneSoup

options(
    DATA_DIR = 'Data',
    OUT_DIR = 'Results',
)

@task
def touch_data():
    for path, _, files in os.walk(options.DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            sh('touch %s' % os.path.join(path, f))


@task
def run():

    ruffus.pipeline_run([top_function])

@task
@needs('touch_data', 'run')
def new_run():
    pass


@ruffus.follows('get_sequences')
def top_function():
    pass

@ruffus.files(os.path.join(options.DATA_DIR, 'ListFiles', 'search_sentinal'),
              os.path.join(options.DATA_DIR, 'ListFiles', 'sequences.list'))
def get_sequence_ids(in_file, out_file):



    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)
    print len(id_list)
    with open(out_file, 'w') as handle:
        handle.write('\n'.join(id_list))


@ruffus.files(os.path.join(options.DATA_DIR, 'ListFiles', 'sequences.list'),
              os.path.join(options.DATA_DIR, 'RawSequences', 'download_sentinal'))
@ruffus.follows(get_sequence_ids)
def get_sequences(in_file, out_file):

    dump_dir = os.path.join(options.DATA_DIR, 'RawSequences')
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

def GetSeqs(client, ID_LIST):

    xml = client.service.run_eFetch(db = 'nucleotide', id = ','.join(ID_LIST))
    soup = BeautifulStoneSoup(xml)
    for seq in soup.findAll('gbseq'):
        gi = None
        for id in seq.findAll('gbseqid'):
            if id.contents[0].startswith('gi'):
                gi = id.contents[0].split('|')[1]
        if gi is None:
            continue

        nt_seq = seq.find('gbseq_sequence').contents[0]
        yield nt_seq, gi




def SearchNCBI(search_sent, recent_date = None, BLOCK_SIZE = 100000, START = 0, db = 'nucleotide'):

    POST_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&'  % db
    POST_URL += 'retmax=%i&' % BLOCK_SIZE
    if START > 0:
        POST_URL += 'retstart=%i&' % START

    search_term = search_sent.replace(' ', '%20')
    search_term = search_term.replace('-', '%20')
    search_term = search_term.replace('+', '%20')
    search_url = POST_URL + '&term=' + search_term

    if recent_date:
        time_delta = datetime.today()-recent_date
        search_url += '&reldate=' + str(time_delta.days)

    xml_data = urllib2.urlopen(search_url).read()

    id_list = re.findall('<Id>(\d*)</Id>', xml_data)
    print len(id_list)
    if len(id_list) >= BLOCK_SIZE-1:
        return id_list + SearchNCBI(search_sent, recent_date = recent_date,
                                      BLOCK_SIZE = BLOCK_SIZE, START = START+BLOCK_SIZE-1)
    else:
        return id_list

def take(N, iterable):
    return list(islice(iterable, N))