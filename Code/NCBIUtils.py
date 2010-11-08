import re, urllib2
from datetime import datetime
from itertools import islice
from BeautifulSoup import BeautifulStoneSoup
from suds.client import Client
import time


def take(N, iterable):
    return list(islice(iterable, N))


def SearchNCBI(search_sent, recent_date = None, BLOCK_SIZE = 1000000, START = 0, db = 'nucleotide'):

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


def extract_sequences(soup, XML = True):

    for seq in soup.findAll('gbseq'):
        gi = None
        for id in seq.findAll('gbseqid'):
            if id.contents[0].startswith('gi'):
                gi = id.contents[0].split('|')[1]
        if gi is None:
            continue

        if XML:
            yield seq.prettify(), gi
        else:
            nt_seq = seq.find('gbseq_sequence').contents[0]
            yield nt_seq, gi


def GetSeqs(ID_LIST, BLOCK_SIZE = 100, XML = False):

    

    wsdl_url = 'http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/efetch_seq.wsdl'
    client = Client('http://www.ncbi.nlm.nih.gov/entrez/eutils/soap/v2.0/efetch_seq.wsdl',
                    retxml = True)
    BAD_GATEWAY = False
    items_left = len(ID_LIST)
    iter_list = iter(ID_LIST)
    block = take(BLOCK_SIZE, iter_list)
    while block:
        try:
            xml = client.service.run_eFetch(db = 'nucleotide', id = ','.join(block))
        except:
            if BAD_GATEWAY:
                break
            time.sleep(60)
            BAD_GATEWAY = True
            continue
        BAD_GATEWAY = False
        soup = BeautifulStoneSoup(xml)
        for seq, gi in extract_sequences(soup, XML):
            yield seq, gi
        items_left -= len(block)
        print 'Items left: %i' % items_left
        block = take(BLOCK_SIZE, iter_list)

def extract_features(in_file):
    
    def get_interval(feature):
        locs = []
        for intervals in feature.findAll('gbfeature_intervals'):
            
            for interval in intervals.findAll('gbinterval'):
                try:
                    locs.append((int(interval.gbinterval_from.contents[0]),
                            int(interval.gbinterval_to.contents[0])))
                except AttributeError:
                    pass
        return locs
        
        
    
    with open(in_file) as handle:
        soup = BeautifulStoneSoup(handle.read())
        
    for feature in soup.findAll('gbfeature'):
        outdict = {'interval':get_interval(feature)}
        
        for qualifier in feature.findAll('gbqualifier'):
            if qualifier.gbqualifier_name.contents[0].strip() == 'product':
                outdict['name'] = qualifier.gbqualifier_value.contents[0].strip()
            elif qualifier.gbqualifier_name.contents[0].strip() == 'translation':
                outdict['AA'] = qualifier.gbqualifier_value.contents[0].strip()
        
        if 'name' in outdict and 'AA' in outdict:
            yield outdict









