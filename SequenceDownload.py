from urllib2 import urlopen
from urlparse import urljoin
from collections import defaultdict
from itertools import izip, repeat, chain, product
import csv, argparse, os.path



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Genomic Sequence Download Code')
    parser.add_argument('--file', dest = 'file', type=str, help='file to check')
    parser.add_argument('--no-download', dest = 'skipdownload', type=bool, action = 'store_true')
    
    args = parser.parse_args()

    base_url = 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria/'
    draft_url = 'ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/'
    wanted_extensions = ('.ptt', '.faa', '.ptt.tgz', '.faa.tgz')
    
    with open(args.file) as handle:
        search_dict = dict(csv.reader(handle))
        
    base_handle = urlopen(base_url)
    draft_handle = urlopen(draft_url)
    iterable = chain(izip(repeat(base_url), csv.reader(base_handle, delimiter = ' ')), 
                    izip(repeat(draft_url), csv.reader(draft_handle, delimiter = ' ')))
    check_locs = defaultdict(list)

    for url, row in iterable:
        if row[0].startswith('d'):        
            for term in search_dict.keys():        
                if row[-1].startswith(term):
                    check_locs[term].append(urljoin(url, row[-1]))

    for term, urls in check_locs.items():
        dump_dir = search_dict[term]
        for url in urls:
            print term, url
            handle = urlopen(url)
            rows = list(csv.reader(handle, delimiter = ' '))
            for row, ext in product(rows, wanted_extensions):
                if row[-1].endswith(ext) and 'scaffold' not in row[-1]:
                    nfile = os.path.join(dump_dir, row[-1])
                    if os.path.exists(nfile):
                        continue
                    nurl = url + '/' + row[-1]
                    ihandle = urlopen(nurl)
                    with open(nfile, 'w') as ohandle:
                        ohandle.write(ihandle.read())
                
    

