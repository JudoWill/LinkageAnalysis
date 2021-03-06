import csv
import os, os.path
from collections import deque
from types import ListType, TupleType, StringType
from itertools import islice, groupby, imap, starmap, repeat, dropwhile, chain
from operator import itemgetter
from functools import partial
from multiprocessing import Lock
from struct import Struct
from fileinput import FileInput


from LinkFields import LINK_FIELDS
import AlignUtils

def AggregateLinkageData(files, full_file, short_file, mode = 'w', short_cut = 0.8):
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
            with open(f) as handle:
                for row in csv.DictReader(handle, delimiter = '\t'):
                    if not row['Total-Num'].startswith('too') and row['Total-Num'] != 'Total-Num':
                        row['Source-Prot'] = sprot
                        row['Target-Prot'] = tprot
                        yield row
    
    files.sort()
    outfields = LINK_FIELDS
    pred = itemgetter('Source-Prot', 'Target-Prot','Source-Start',
                     'Source-End', 'Target-Start','Target-End')
    
    if mode == 'w':
        item_iter = multi_file_iterator(files)
    elif mode == 'a':
        try:
            with open(full_file) as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                items = deque(reader, 1)
                last_item = items.pop()
                lvals = pred(last_item)
                item_iter = dropwhile(lambda x: pred(x) <= lvals, 
                                            multi_file_iterator(files))


        except IOError:
            item_iter = multi_file_iterator(files)
        except IndexError:
            item_iter = multi_file_iterator(files)
            
    
    total_getter = itemgetter('Total-Num')
    
    with open(full_file, mode) as fhandle:
        fwriter = csv.DictWriter(fhandle, outfields, delimiter = '\t', extrasaction = 'ignore')
        if mode == 'w':
            fwriter.writerow(dict(zip(outfields, outfields)))
        with open(short_file, mode) as shandle:
            swriter = csv.DictWriter(shandle, outfields, delimiter = '\t', extrasaction = 'ignore')
            if mode == 'w':
                swriter.writerow(dict(zip(outfields, outfields)))
            
            for key, group in groupby(item_iter, pred):
                lgroup = list(group)

                
                if float(lgroup[0]['Total-Score']) >= short_cut:
                    fwriter.writerows(lgroup)
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
    a1_mapping, _ = aln1.convert_numbering(ref_genome)
    #print len(_), len(a1_numbering)
    #raise KeyError
    aln2 = AlignUtils.Alignment.alignment_from_file(afile2)
    a2_mapping, _ = aln2.convert_numbering(ref_genome)

    #a1_mapping = dict(zip(range(len(a1_numbering)), a1_numbering))
    #a2_mapping = dict(zip(range(len(a2_numbering)), a2_numbering))

    #fall back numbers when we fall out of range
    fb1 = max(a1_mapping)
    fb2 = max(a2_mapping)

    conv_fields = (('Source-Start', a1_mapping, max(a1_mapping)),
                    ('Source-End', a1_mapping, max(a1_mapping)),
                    ('Target-Start', a2_mapping, max(a2_mapping)), 
                    ('Target-End', a2_mapping, max(a2_mapping)))


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
                        row[field] = mapping[int(row[field])]
                    except ValueError:
                        skip = True
                        break
                    except IndexError:
                        print 'too large', int(row[field]), len(mapping)
                        row[field] = fb
                if not skip:
                    writer.writerow(row)

def merge_to_binary(indirec, outfile):
    

    def safefloat(string):
        try:
            return float(string)
        except ValueError:
            if len(string) == 0:
                return 0.0
            else:
                return None

    def safeint(string):
        try:
            return int(string.strip())
        except ValueError:
            if len(string) == 0:
                return 0
            else:
                return None


    files = sorted([x for x in os.listdir(indirec) if x.endswith('.conv')])
    max_len = max(len(x.split('--')[0]) for x in files)

    mapping = (('Start','i', safeint), ('End','i', safeint), ('Num', 'i', safeint),
                ('Prot', '%ic'%max_len, lambda x: x.ljust(max_len)),
                ('Score', 'f', safefloat), ('Cons', 'f', safefloat), 
                ('Info', 'f', safefloat), ('Dist', 'f', safefloat),
                ('Seq', 'c', lambda x:x))

    mapping_list = []
    fmtstr = ''
    for field in LINK_FIELDS:
        for key, fmt, func in mapping:
            if key in field:
                print field, key, fmt
                fmtstr += fmt
                mapping_list.append(func)
                break
    StructClass = Struct(fmtstr)

    
    handle = FileInput([os.path.join(indirec, x) for x in files])
    reader = csv.reader(handle, delimiter = '\t')
    grouper = itemgetter(*range(6)) #group by source-prot through target-end
    buf = open(outfile, 'wb')
    count = 0
    for ind, (key, rows) in enumerate(groupby(reader, key = grouper)):
        trow = rows.next()
        converted = [func(field) for func, field in zip(mapping_list, trow)]
        #print zip(LINK_FIELDS, converted, trow)
        if all(x is not None for x in converted):
            arglist = list(chain(converted[0], converted[1], converted[2:]))
            data = StructClass.pack(*arglist)
            buf.write(data)
            count += 1
            if count % 10000 == 0:
                print count, key


def greatest_line(filename):
    
    with open(filename) as handle:
        greatest = None
        for line in handle:
            if line > greatest:
                greatest = line
    return greatest


class LockedFile(file):
    def __init__(self,name,mode):
        self.lock = Lock()
        file.__init__(self,name,mode)

    def safe_write(self, *args, **kwargs):
        with self.lock:
            self.write(*args, **kwargs)


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




