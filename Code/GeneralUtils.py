import csv
import os, os.path
from collections import deque
from types import ListType, TupleType
from itertools import islice, groupby, imap, starmap, repeat
from operator import itemgetter

def unique_justseen(iterable, key = None):
    """Yields unique values."""
    return imap(next, imap(itemgetter(1), groupby(iterable, key)))

def repeatfunc(func, times , *args):
    """Repeatedly causes a function."""
    if times is None:
        return starmap(func, repeat(args))
    return starmap(func, repeat(args, times))

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
    """Takes N items from an iterable."""
    return list(islice(iterable, N))
    
    
def OverlappingIterator(iterable, win_size, win_overlap):
    """Yields overlapping sets of items."""
    
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
    """Iterates over a fasta-file

    Yields (name, seq) pairs."""

    name = None
    with open(filename) as handle:
        for header, group in groupby(handle, lambda x: x.startswith('>')):
            if header:
                name = group.next().strip()[1:]
            else:
                seq = ''.join([x.strip() for x in group])
                yield name, seq

def count_fasta(filename):
    """Counts the number of sequences in a fasta-file."""

    return sum(imap(bool, fasta_iter(filename)))
    
def join_fasta(filenames, out_file, mode = 'w', strip = False):
    """Joins multiple fasta file into one."""
    
    with open(out_file, mode) as ohandle:
        for f in filenames:
            for name, seq in fasta_iter(f):
                if strip:
                    name = gi_from_path(name)
                    name = name.split('_')[0]
                ohandle.write('>%s\n%s\n' % (name, seq))
            
def split_fasta(filename, max_num = 20000):
    """Splits fasta files."""
    
    it = fasta_iter(filename)
    total_num = sum(imap(bool, it))
    print total_num
    if total_num < max_num:
        return None
    
    out_files = []
    it = fasta_iter(filename)
    num = 0
    nseqs = max_num-2
    while nseqs+2 >= max_num:
        num += 1
        fname = filename+str(num)
        out_files.append(fname)
        with open(fname, 'w') as handle:
            for nseqs, (name, seq) in enumerate(islice(it, max_num-1)):
                handle.write('>%s\n%s\n' % (name, seq))

    return out_files



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
    """Returns the GI from a path."""

    fname = path.split(os.sep)[-1]
    gi = fname.split('.')[0]
    return gi
    
def prots_from_path(path):
    
    fname = path.split(os.sep)[-1]
    p1, p2 = fname.split('.')[0].split('--')
    return p1,p2
    
class pushd():
    """Chnages a directory and then changes back when leaving the context."""

    def __init__(self, newpath):
        self.prev_path = os.getcwd()
        self.new_path = newpath        

    def __enter__(self):
        os.chdir(self.new_path)
    def __exit__(self, typ, value, tb):
        os.chdir(self.prev_path)        



def safe_mkdir(path):
    """Makes a new directory but catches Exceptions."""
    
    try:
        os.makedirs(path)
    except OSError:
        pass




