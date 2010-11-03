from paver.easy import *
import csv, os, os.path
import ruffus
import urllib2, re
from datetime import datetime

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


@ruffus.follows('get_sequence_ids')
def top_function():
    pass

@ruffus.files(os.path.join(options.DATA_DIR, 'ListFiles', 'search_sentinal'),
              os.path.join(options.DATA_DIR, 'ListFiles', 'sequences.list'))
def get_sequence_ids(in_file, out_file):



    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)

    with open(out_file, 'w') as handle:
        handle.write('\n'.join(id_list))


@ruffus.files(os.path.join(options.DATA_DIR, 'ListFiles', 'sequences.list'),
              os.path.join(options.DATA_DIR, 'RawSequences', 'download_sentinal'))
def get_sequences(in_file, out_file):
    pass







######Utility Functions!

def GetXML(ID_LIST, db = 'nucleotide'):

    POST_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db=%s' % db
    RET_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=%s&query_key=1&mode=xml&rettype=full' % db

    pmid_list = ','.join(map(lambda x: str(x), ID_LIST))
    post_req_url = POST_URL + '&id=' + pmid_list

    post_res = urllib2.urlopen(post_req_url).read()

    web_env = re.findall('<WebEnv>(.*?)</WebEnv>', post_res)[0]

    req_url = RET_URL + '&WebENV=' + web_env

    xml_data = urllib2.urlopen(req_url).read()
    return xml_data.decode('ascii', 'ignore')


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

    if len(id_list) == BLOCK_SIZE:
        return id_list + SearchNCBI(search_sent, recent_date = recent_date,
                                      BLOCK_SIZE = BLOCK_SIZE, START = START+BLOCK_SIZE-1)
    else:
        return id_list