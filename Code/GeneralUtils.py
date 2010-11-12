import csv
import os
from collections import deque
from types import ListType, TupleType
from itertools import islice

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
