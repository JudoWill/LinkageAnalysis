from urllib2 import urlopen
from urlparse import urljoin
from collections import defaultdict
from itertools import izip, repeat, chain, product, groupby
from operator import itemgetter
from functools import partial
import csv, argparse, os.path, logging, tarfile, os
from shutil import rmtree


def get_directories(url, check_terms):
    """Yields the directories from an FTP-url."""

    handle = urlopen(url)
    for row in csv.reader(handle, delimiter = ' '):
        for term in check_terms:
            if row[0].startswith('d') and row[-1].startswith(term):
                yield term, url + '/' + row[-1]

def get_file_urls(url):
    """Yields the file urls that should be downloaded."""

    wanted_extensions = ('.ptt', '.faa', '.ptt.tgz', '.faa.tgz')
    handle = urlopen(url)
    filenames = map(itemgetter(-1), csv.reader(handle, delimiter = ' '))
    for fname, ext in product(filenames, wanted_extensions):
        if fname.endswith(ext) and 'scaffold' not in fname:
            yield url + '/' + fname

def download_file(url, dump_dir):
    """Downloads a url and puts into teh dump_dir"""
    
    fname = url.rsplit('/',1)[-1]
    with open(os.path.join(dump_dir, fname), 'w') as handle:
        handle.write(urlopen(url).read())

def check_NCBI(search_dict, urls):
    """Checks the NCBI ftp directories and downloads the file."""

    check_locs = defaultdict(list)
    logging.warning('Getting base directories')
    for base_url in urls:
        for term, url in get_directories(base_url, search_dict.keys()):
            check_locs[term].append(url)
    
    for term, urls in check_locs.items():
        dump_dir = search_dict[term]
        try:
            os.makedirs(dump_dir)
        except OSError:
            pass

        logging.warning('Fetching %i genomes for %s' % (len(urls), term))
        for url in urls:
            logging.warning('Fetching: ' + url)
            for furl in get_file_urls(url):
                download_file(furl, dump_dir)

def unzip_dir(base_dir):
    """Unzips all of the relevant files."""
    
    flist = [x for x in os.listdir(base_dir) if x.endswith('.tgz')]
    for f in flist:
        fname = os.path.join(base_dir, f)
        handle = tarfile.open(fname)
        handle.extractall(path = base_dir)
        handle.close()

def read_ptt(filename):
    """Reads ptt files and returns a name-to-gi mapping."""
    
    with open(filename) as handle:
        handle.next()
        handle.next()
        id2name = {}
        for row in csv.DictReader(handle, delimiter = '\t'):
            if row['Gene'] != '-':            
                id2name[row['PID']] = row['Gene']
    return id2name            

def read_fasta(filename):
    """Reads a fasta file and returns a go-to-seq mapping."""
    
    id2seq = {}
    with open(filename) as handle:
        name = None
        for key, lines in groupby(handle, lambda x: x.startswith('>')):
            if key:
                name = lines.next().strip()[1:]
                gi = name.split('|')[1]
            else:
                id2seq[gi] = ''.join([x.strip() for x in lines])
    return id2seq

def process_pair(fasta_file, ptt_file):
    """Process .faa and .ptt files to return name-to-seq mapping."""
    
    gi2seq = read_fasta(fasta_file)
    gi2name = read_ptt(ptt_file)

    name2seq = dict()
    for gi, seq in gi2seq.iteritems():
        if gi in gi2name:        
            name2seq[gi2name[gi]] = seq

    return name2seq

def process_many(grouped_names):
    """Processes groups of .faa and .ptt files."""

    res = dict()    
    for group in grouped_names:
        try:        
            res.update(process_pair(group + '.faa', group + '.ptt'))
        except IOError:
            pass

    return res

def process_directory(direc):
    
    files = os.listdir(direc)
    df = partial(os.path.join, direc)
    files = [df(x.split('.')[0]) for x in files if x.endswith('.faa')]
    grouper = lambda x: x.split(os.sep)[-1][:9]
    results = []
    logging.warning('%i files to process' % len(files))
    for key, fs in groupby(sorted(files), grouper):
        logging.warning('Processing group: ' + key)
        results.append((key, process_many(fs)))
    
    counts = defaultdict(int)
    for _, res in results:
        for key in res.iterkeys():
            counts[key] += 1

    for key, count in counts.iteritems():
        if count >= 10:
            try:
                fname = df('Aggregated', key + '.fasta')
                with open(fname, 'w') as handle:
                    for genome, seq_dict in results:
                        if key in seq_dict:
                            handle.write('>%s\n%s\n' % (genome, seq_dict[key]))
            except IOError:
                pass


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Genomic Sequence Download Code')
    parser.add_argument('--file', dest = 'file', help='file to check')
    parser.add_argument('--no-download', dest = 'skipdownload', default = False,
                        action = 'store_true')
    parser.add_argument('--no-drafts', dest = 'skipdraft', default = False,
                        action = 'store_true')
    parser.add_argument('--no-unzip', dest = 'skipunzip', default = False,
                        action = 'store_true')
    parser.add_argument('--fresh', dest = 'fresh', default = False,
                        action = 'store_true')
    parser.add_argument('--min-required', dest = 'minrequired', default = 7)
    
    args = parser.parse_args()
    
    urls = ('ftp://ftp.ncbi.nih.gov/genomes/Bacteria/',)
    if not args.skipdraft:
        urls += ('ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/',)
    
    with open(args.file) as handle:
        search_dict = dict(csv.reader(handle))

    if args.fresh:
        logging.warning('Cleaning up previous downlaods!')
        for d in search_dict.values():
            rmtree(d)
            os.mkdir(d)
        
        
    if not args.skipdownload:
        check_NCBI(search_dict, urls)
    
    if not args.skipunzip:
        for direc in search_dict.values():
            logging.warning('Unzipping: ' + direc)
            unzip_dir(direc)
    
    for key, direc in search_dict.items():
        logging.warning('Processing ' + key)
        try:
            os.makedirs(os.path.join(direc, 'Aggregated'))
        except OSError:
            pass
        process_directory(direc)
















