import csv
import os, os.path
from collections import deque
from types import ListType, TupleType, StringType
from itertools import islice, groupby, imap, starmap, repeat, dropwhile
from operator import itemgetter
from functools import partial

import AlignUtils

def AggregateLinkageData(files, full_file, short_file, mode = 'w', short_cut = 0.9):
    """Aggregates linkage data into to files for easier processing.
    
    Inputs:
    in_direc    The input directory which has all linkages.
    full_file   The place to put all aggregaed data.
    short_file  The place to put all shortened data.
    """
    
    def multi_file_iterator(files):
        for f in files:
            fpart = f.split(os.sep)[-1]
            sprot, tprot = fpart.split('.')[0].split('--')
            with open(os.path.join(in_direc,f)) as handle:
                for row in csv.DictReader(handle, delimiter = '\t'):
                    if not row['Total-Num'].startswith('too') and row['Total-Num'] != 'Total-Num':
                        row['Source-Prot'] = sprot
                        row['Target-Prot'] = tprot
                        yield row
    
    files.sort()
    outfields = ('Source-Prot', 'Target-Prot','Source-Start','Source-End',
                'Target-Start','Target-End','Source-Seq','Target-Seq',
                'Correct-Num','Total-Num','This-Score','Total-Score')
    pred = itemgetter('Source-Prot', 'Target-Prot','Source-Start',
                     'Source-End', 'Target-Start','Target-End')
    
    if mode == 'w':
        item_iter = multi_file_iterator(files)
    elif mode == 'a':
        with open(full_file) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            items = deque(reader, 1)
            try:
                last_item = items.pop()
                lvals = pred(last_item)
                item_iter = dropwhile(lambda x: pred(x) <= lvals, 
                                        multi_file_iterator(files, in_direc))
            except IndexError:
                item_iter = multi_file_iterator(files)
    
    total_getter = itemgetter('Total-Num')
    
    with open(full_file, mode) as fhandle:
        fwriter = csv.DictWriter(fhandle, outfields, delimiter = '\t')
        if mode == 'w':
            fwriter.writerow(dict(zip(outfields, outfields)))
        with open(short_file, mode) as shandle:
            swriter = csv.DictWriter(shandle, outfields, delimiter = '\t')
            if mode == 'w':
                swriter.writerow(dict(zip(outfields, outfields)))
            
            for key, group in groupby(item_iter, pred):
                lgroup = list(group)
                fwriter.writerows(lgroup)
                if float(lgroup[0]['Total-Score']) >= short_cut:
                    total_seqs = sum((int(x['Total-Num']) for x in lgroup))
                    total_correct = sum((int(x['Correct-Num']) for x in lgroup))
                    lgroup[0]['Total-Num'] = total_seqs
                    lgroup[0]['Correct-Num'] = total_correct
                    swriter.writerow(lgroup[0])
                    sync_handle(shandle)
                sync_handle(fhandle)
                

def convert_numbering(afile1, afile2, link_file, out_file, ref_genome):
    """Converts the numbering in a linkage file to the numbering in the reference genome."""

    aln1 = AlignUtils.Alignment.alignment_from_file(afile1)
    _, a1_numbering = aln1.convert_numbering(ref_genome)
    aln2 = AlignUtils.Alignment.alignment_from_file(afile2)
    _, a2_numbering = aln2.convert_numbering(ref_genome)

    a1_mapping = dict(zip(range(len(a1_numbering)), a1_numbering))
    a2_mapping = dict(zip(range(len(a2_numbering)), a2_numbering))

    #fall back numbers when we fall out of range
    fb1 = len(a1_numbering)
    fb2 = len(a2_numbering)

    conv_fields = (('Source-Start', a1_mapping, len(a1_numbering)),
                    ('Source-End', a1_mapping, len(a1_numbering)),
                    ('Target-Start', a2_mapping, len(a2_numbering)), 
                    ('Target-End', a2_mapping, len(a2_numbering)))


    with open(link_file) as handle:
        fields = handle.next().split('\t')

    with open(link_file) as handle:
        reader = csv.DictReader(handle, delimiter = '\t')
        with open(out_file, 'w') as handle:
            writer = csv.DictWriter(handle, fields, delimiter = '\t', extrasaction = 'ignore')
            writer.writerow(dict(zip(fields, fields)))
            for row in reader:
                skip = False
                for field, mapping, fb in conv_fields:
                    try:
                        row[field] = mapping.get(int(row[field]), fb)
                    except ValueError:
                        skip = True
                        break
                if not skip:
                    writer.writerow(row)



def touch_existing(fnames, times = None):
    """Touches existing files."""

    if type(fnames) == StringType:
        fnames = [fnames]

    for fname in fnames:
        if os.path.exists(fname):
            with file(fname, 'a'):
                os.utime(fname, times)


            
def sync_handle(fhandle):
    """Helper function for syncing filehandles to disk."""
    fhandle.flush()
    os.fsync(fhandle.fileno())
    
    
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




