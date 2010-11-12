import csv
import os, os.path
from collections import deque
from types import ListType, TupleType
from itertools import islice, groupby, imap

def filter_gi(load_dir, dump_dir, load_extension = '.xml', 
                dump_extension = '.xml', force_new = False):
    
    done = set()
    if not force_new:
        for f in os.listdir(dump_dir):
            if f.endswith(load_extension):
                done.add(gi_from_path(f))
    need = set()
    for f in os.listdir(load_dir):
        if f.endswith(load_extension):
            need.add(gi_from_path(f))
    
    items = need-done
    for count, gi in enumerate(items):
        if count % 500 == 0:
            print 'Processing: %i of %i' % (count, len(items))
        yield gi


def take(N, iterable):
    return list(islice(iterable, N))
    
    
def OverlappingIterator(iterable, win_size, win_overlap):
    
    this_iter = iter(iterable)
    
    items = deque([], win_size)
    items.extend(take(win_size, this_iter))
    
    yield items
    
    block = take(win_overlap, this_iter)
    while block:
        items.extend(block)
        yield items
        block = take(win_overlap, this_iter)
        
    
    
def fasta_iter(filename):
    name = None
    with open(filename) as handle:
        for header, group in groupby(handle, lambda x: x.startswith('>')):
            if header:
                name = group.next().strip()
            else:
                seq = ''.join([x.strip() for x in group])
                yield name, seq
    
    
    

def make_mapping_dict(in_file):
    
    mapping_dict = {}
    with open(in_file) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            if row['name'] == 'None':
                mapping_dict[row['key']] = None
            else:
                mapping_dict[row['key']] = row['name']
    return mapping_dict


def mapping_func(mapping_dict, name):
    return mapping_dict.get(name.lower(), None)

def gi_from_path(path):

    fname = path.split(os.sep)[-1]
    gi = fname.split('.')[0]
    return gi
