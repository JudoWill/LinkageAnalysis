import re, urllib2
from datetime import datetime
from itertools import islice
from BeautifulSoup import BeautifulStoneSoup, SoupStrainer
from suds.client import Client
import time
from collections import defaultdict
from subprocess import call
import shlex
import re
from xml.etree.ElementTree import ElementTree
from xml.parsers.expat import ExpatError

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


def extract_sequences(soup, XML = True, seq_only = False):

    for seq in soup.findAll('gbseq'):
        gi = None
        if not seq_only:        
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

def extract_features(in_file, mapping = None):
    
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
                name = qualifier.gbqualifier_value.contents[0].strip().lower()
                if mapping is not None:
                    name = mapping(name)
                if name is None:
                    break
                outdict['name'] = name
            elif qualifier.gbqualifier_name.contents[0].strip() == 'translation':
                outdict['AA'] = qualifier.gbqualifier_value.contents[0].strip()
        
        if 'name' in outdict and 'AA' in outdict:
            yield outdict


def determine_subtype(in_file):
    hits = defaultdict(int)
    with open(in_file) as handle:
        soup = BeautifulStoneSoup(handle.read())
    
    for seq in soup.findAll('iteration'):
        try:
            hit = seq.iteration_hits.hit.hit_def.contents[0]
        except:
            hit = None
        if hit:
            hits[hit.split('_')[1]] += 1
    
    count = sum(hits.values())
    if count < 5:
        return None
    elif all([x < count*0.6 for x in hits.values()]):
        #print 'too heterogenus %s' % ','.join(map(str,hits.items()))
        return None
    else:
        for key, val in hits.items():
            if val > count*0.6:
                return key

def determine_subtype_short(in_file):
    hits = defaultdict(int)
    strainer = SoupStrainer(re.compile('iteration'))
    with open(in_file) as handle:
        soup = BeautifulStoneSoup(handle.read(), parseOnlyThese = strainer)
    
    for seq in soup.findAll('iteration'):
        try:
            hit = seq.iteration_hits.hit.hit_def.contents[0]
        except:
            hit = None
        if hit:
            hits[hit.split('_')[1]] += 1
    
    count = sum(hits.values())
    if count < 5:
        return None
    elif all([x < count*0.6 for x in hits.values()]):
        print 'too heterogenus %s' % ','.join(map(str,hits.items()))
        return None
    else:
        for key, val in hits.items():
            if val > count*0.6:
                return key

def determine_subtype_element(in_file):
    hits = defaultdict(int)
    try:    
        tree = ElementTree(file = in_file)

        for it in tree.getiterator('Iteration'):
            hit_list = it.getiterator('Hit')
            if len(hit_list) > 0:
                hit = hit_list[0].find('Hit_def').text
                hits[hit.split('_')[1]] += 1
    except ExpatError:
        return None
    
    count = sum(hits.values())
    if count < 5:
        return None
    elif all([x < count*0.6 for x in hits.values()]):
        print 'too heterogenus %s' % ','.join(map(str,hits.items()))
        return None
    else:
        for key, val in hits.items():
            if val > count*0.6:
                return key
        
            
def make_blast_cmd(program_type, database_path, in_path, out_path, blast_type = 0, **options):
    
    def get_formatdb_options(options):
        assert 'dbtype' in options, \
                'You must have a "dbtype" option when using formatdb'
        dbtype = options['dbtype'].lower()
        assert dbtype.startswith('p') or dbtype.startswith('n'), \
            'Arguement to dbtype must either be "prot" or "nuc"'
        return dbtype
        
        
    
    if blast_type == 0:
        if program_type == 'formatdb':
            dbtype = get_formatdb_options(options)
            if dbtype.startswith('p'):
                prot = 'T'
            elif dbtype.startswith('n'):
                prot = 'F'
            info = {
            'ipath':in_path,
            'prot':prot
            } 
            return 'formatdb -i %(ipath)s -p %(prot)s' % info
        elif program_type == 'blastn':
            info = {
            'dpath':database_path,
            'ipath':in_path,
            'opath':out_path
            }
            return 'blastall -p blastn -d %(dpath)s -i %(ipath)s -m 7 -o %(opath)s' % info
    elif blast_type == 1:
        if program_type == 'formatdb':
            dbtype = get_formatdb_options(options)
            if dbtype.startswith('p'):
                prot = 'prot'
            elif dbtype.startswith('n'):
                prot = 'nucl'
            info = {
            'ipath':in_path,
            'prot':prot,
            }
            return 'makeblastdb -in %(ipath)s -dbtype %(prot)s' % info
        elif program_type == 'blastn':
            info = {
            'dpath':database_path,
            'ipath':in_path,
            'opath':out_path
            }
            return 'blastn -db %(dpath)s -query %(ipath)s -out %(opath)s -outfmt 5' % info
        elif program_type = 'blastx':
            info = {
            'dpath':database_path,
            'ipath':in_path,
            'opath':out_path
            }
            return 'blastx -db %(dpath)s -query %(ipath)s -out %(opath)s -outfmt 5 -best_hit_overhang 0.25' % info
            

def guess_blast_computer_type():
    
    progs = ('formatdb', 'makeblastdb')
    for i, prog in enumerate(progs):
        print i
        args = shlex.split('which ' + prog)        
        retcode = call(args)
        if retcode == 0:
            return i        
    assert False, 'Could not find BLAST on this computer!'



