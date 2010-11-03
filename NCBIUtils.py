import re, urllib2
from datetime import datetime
from itertools import islice


def take(N, iterable):
    return list(islice(iterable, N))


def SearchNCBI(search_sent, recent_date = None, BLOCK_SIZE = 100000, START = 0, db = 'nucleotide'):

    POST_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=%s&'
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

    if len(id_nums) == BLOCK_SIZE:
        return id_list + SearchNCBI(search_sent, recent_date = recent_date,
                                      BLOCK_SIZE = BLOCK_SIZE, START = START+BLOCK_SIZE-1)
    else:
        return id_list
