from urllib2 import urlopen
from urlparse import urljoin
from collections import defaultdict
from itertools import izip, repeat, chain, product
from operator import itemgetter
import csv, argparse, os.path, logging, tarfile



def get_directories(url, check_terms):
    handle = urlopen(url)
    for row in csv.reader(handle, delimiter = ' '):
        for term in check_terms:
            if row[0].startswith('d') and row[-1].startswith(term):
                yield term, url + '/' + row[-1]

def get_file_urls(url):
    wanted_extensions = ('.ptt', '.faa', '.ptt.tgz', '.faa.tgz')
    handle = urlopen(url)
    filenames = map(itemgetter(-1), csv.reader(handle, delimiter = ' '))
    for fname, ext in product(filenames, wanted_extensions):
        if fname.endswith(ext) and 'scaffold' not in fname:
            yield url + '/' + fname

def download_file(url, dump_dir):
    
    fname = url.rsplit('/',1)[-1]
    with open(os.path.join(dump_dir, fname), 'w') as handle:
        handle.write(urlopen(url).read())

def check_NCBI(search_dict, urls):

    check_locs = defaultdict(list)
    logging.warning('Getting base directories')
    for base_url in urls:
        for term, url in get_directories(base_url, search_dict.keys()):
            check_locs[term].append(url)
    
    for term, urls in check_locs.items():
        dump_dir = search_dict[term]
        logging.warning('Fetching %i genomes for %s' % (len(urls), term))
        for url in urls:
            logging.warning('Fetching: ' + url)
            for furl in get_file_urls(url):
                download_file(furl, dump_dir)

def unzip_dir(base_dir):
    
    flist = [x for x in os.listdir(base_dir) if x.endswith('.tgz')]
    for f in flist:
        fname = os.path.join(base_dir, f)
        handle = tarfile.open(fname)
        handle.extractall()
        handle.close()
        


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Genomic Sequence Download Code')
    parser.add_argument('--file', dest = 'file', help='file to check')
    parser.add_argument('--no-download', dest = 'skipdownload', default = False,
                        action = 'store_true')
    parser.add_argument('--no-drafts', dest = 'skipdraft', default = False,
                        action = 'store_true')
    parser.add_argument('--no-unzip', dest = 'skipunzip', default = False,
                        action = 'store_true')
    
    args = parser.parse_args()

    urls = ('ftp://ftp.ncbi.nih.gov/genomes/Bacteria/',)
    if not args.skipdraft:
        urls += ('ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/',)
    
    with open(args.file) as handle:
        search_dict = dict(csv.reader(handle))
        
    if not args.skipdownload:
        check_NCBI(search_dict, urls)
    
    if not args.skipunzip:
        for direc in search_dict.values():
            logging.warning('Unzipping: ' + direc)
            unzip_dir(direc)
    
















