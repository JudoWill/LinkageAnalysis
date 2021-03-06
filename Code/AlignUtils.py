from __future__ import division
from GeneralUtils import *
from subprocess import call
import shlex
from itertools import groupby, product, dropwhile, imap, izip, count
from math import log
from random import  sample, randint
from operator import itemgetter, ne, eq
import csv
from functools import partial
import  os
from collections import defaultdict
from LinkFields import LINK_FIELDS
from LinkUtils import calculate_mutual_info, calculate_PNAS, prediction_mapping
from LinkUtils import calculate_OMES, calculate_SBASC, get_sub_mat
from LinkUtils import LinkCalculator
from pylru import lrudecorator

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

def take(N, iterable):
    """Takes N items from an iterable."""
    return list(islice(iterable, N))

@lrudecorator(1000)
def calculate_entropy(seq):
    """Salculates the Shannon Entropy of any sequence.
    
    Arguements:
    seq -- Any interable sequence.
    
    Returns:
    float -- entropy of the sequence.
    """

    counts = defaultdict(int)
    for v in seq:
        counts[v] += 1

    ent = 0.0    
    for val in counts.values():
        ent += (val/len(seq))*log(val/len(seq))

    return ent

def make_counts(signal):
    cdict = defaultdict(int)
    for s in signal:
        cdict[s] += 1
    return cdict

class Alignment():
    """A container for holding multiple alignments."""
    
    def __init__(self):
        self.seqs = {}
        self.width = None
        self.seq_nums = []
    
    @staticmethod
    def alignment_from_file(filename):
        """Loads alignments from a tab-delimited file."""

        align = Alignment()        
        with open(filename) as handle:
            for line in handle:
                parts = line.strip().split('\t')
                align.seqs[parts[0]] = parts[1]
                align.width = len(parts[1])
        align._process_seq_nums()
        return align

    def _process_seq_nums(self):
        """Calculating the number of letters at a position."""
        
        self.seq_nums = []
        for group in zip(*self.seqs.values()):
            self.seq_nums.append(len(set(group)))

    def append_alignment(self, newaln):
        """Appends a new alignment onto this one. Removes all keys not present in both."""
        
        valid_keys = set(self.seqs.keys()) & set(newaln.seqs.keys())

        rmkeys = set(self.seqs.keys()) - valid_keys
        for key in rmkeys:
            self.seqs.pop(key)

        for key in valid_keys:
            self.seqs[key] += newaln.seqs[key]

        self.width = len(self.seqs.values()[0])
        self._process_seq_nums()

    def entropy_filter(self, topN, allowed_gaps = 0):
        """Filters out all but the topN entropy sequences."""

        entropies = []
        for group, ind in izip(zip(*self.seqs.values()), count(0)):
            if sum([x == '-' for x in group]) <= allowed_gaps:
                entropies.append((ind, calculate_entropy(group)))
        
        scols = sorted(entropies, key = itemgetter(1))
        wanted_cols = sorted([x for x, _ in scols[:topN]])
        getter = itemgetter(*wanted_cols)
        
        for key, seq in self.seqs.items():
            self.seqs[key] = ''.join(getter(seq))
            self.width = len(self.seqs[key])
        self._process_seq_nums()

    def bootstrap_columns(self, num_reps):
        """Yields randomized bootstrap alignments"""

        for n in xrange(num_reps):
            inds = [randint(0, self.width-1) for _ in xrange(self.width)]
            getter = itemgetter(*inds)
            nalign = Alignment()
            for key, seq in self.seqs.iteritems():
                nalign.seqs[key] = ''.join(getter(seq))
            nalign.width = self.width
            nalign._process_seq_nums()
            yield nalign
            


    def write_phylip(self, fname):
        """Writes a phylip formatted alignment file"""

        with open(fname, 'w') as handle:
            handle.write('%i %i\n' % (len(self.seqs), self.width))

            for key, seq in self.seqs.items():
                handle.write('%s%s\n' % (key[:10].ljust(10), seq))

    def write_fasta(self, fname):
        """Writes the alignment in fasta format"""

        with open(fname, 'w') as handle:
            for key, seq in self.seqs.items():
                iterable = iter(seq)
                handle.write('>%s\n' % key)
                block = take(80, iterable)
                while block:
                    handle.write(''.join(block) + '\n')
                    block = take(80, iterable)

    def write_aln(self, fname):
        """Writes the alignment in ALN format"""

        with open(fname, 'w') as handle:
            for key, seq in self.seqs.items():
                handle.write('%s\t%s\n' % (key, seq))            

    def iterate_columns(self, seq_names):

        tseqs = [self.seqs[key] for key in seq_names]
        for group in izip(*tseqs):
            yield ''.join(group)


    def get_slice(self, start, stop, MIN_NUM = 2):
        """Return a slice of an alignment.

            Returns a sub-slice of an alignment as an Alignment instance. Will 
            return None is the alignment is too conserved.
            
            Arguments:
            start -- Start of the region.
            stop -- End of the region.
            
            Kwargs:
            MIN_NUM -- The minumum number of unique sequences.

            Returns:
            An Alignemnt-instance
            """
        
        if any([x >= MIN_NUM for x in self.seq_nums[start:stop]]):
            ngapfun = partial(ne, '-')
            gapfun = partial(eq, '-')
            align = Alignment()
            for name, seq in self.seqs.items():
                if any(imap(ngapfun, seq[:stop])) and not all(imap(gapfun, seq[start:stop])):
                    align.seqs[name] = seq[start:stop]
            return align
        else:
            return None

    def get_signal(self, seq_names):
        """Returns an interger "signal" from the alignment.

        Arguements:
        seq_names -- An ordered list of sequences that are part of this signal.

        Returns:
        (signal_list, mapping_dict)
        signal_list -- List of integers indicating the signal.
        mapping_dict -- A dictionary which maps the sequence keys 
                        to signal values.
        """
        
        signal = []
        seq_dict = {}
        for name in seq_names:
            seq = self.seqs[name]
            signal.append(seq_dict.setdefault(seq, len(seq_dict)+1))

        return signal, seq_dict

    def get_consensus(self):
        """Returns the consensus sequence.

        Returns a sequence which is the most common value at each position. 
        Gaps are excluded from the count.
        """

        seqs = []
        for col in range(self.width):
            cdict = defaultdict(int)            
            for s in self.seqs.values():
                cdict[s[col]] += 1
            #cdict.pop('-')
            m_item = max(cdict.items(), key = itemgetter(1))
            if m_item[1] > 0.6*sum(cdict.values()) and sum(cdict.values()) > 0.6*len(self.seqs):
                seqs.append(m_item[0])
        return ''.join(seqs)
    
    def convert_numbering(self, dest_key):
        """Determines the numbering between an alignment and a sequence.

        Determines the coulmbn number mapping between columns in an 
        alignment and sequence positions. This allows one to convert results 
        that are in "alignment-space" to a "reference sequence space".
        
        Arguements:
        dest_key -- The desired reference sequence. MUST be present in the 
                    alignment.
        
        Returns:
        (dest_nums, align_nums)
        dest_nums -- A list which is the same length as the "reference" 
                    sequence. The items indicate the corresponding column in 
                    the alignment.
        align_nums -- A list the same length as the alignment. The items 
                      indicate the corresponding position in the reference 
                      sequence."""
        
        seq_count = -1
        align_nums = [] #same length as sequence and tells which column the reference belongs too
        dest_nums = [] #same length as alignment and tells the position in the reference
        for num, let in enumerate(self.seqs[dest_key]):
            if let != '-':
                seq_count += 1
                align_nums.append(num)
            dest_nums.append(seq_count)
            
        return dest_nums, align_nums


def prot_from_path(path):
    """Returns the GI from a path."""

    fname = path.split(os.sep)[-1]
    gi = fname.split('.')[0]
    return gi




def run_muscle(filename, out_align, MAX_MEM = 1500):
    """Runs MUSCLE command-line alignment program.

    Uses subprocessing to call the MUSCLE alignment tool to align sequences 
    in a fasta-file.

    Arguements:
    filename -- A full-path to a fasta-formatted sequence file.
    out_align -- A full-path to the desired dest of the output alignment.

    Kwargs:
    MAX_MEM -- The maximum amount of memory allowed for the MUSCLE program
                in mb. Default: 1500

    Returns:
    None"""    

    tfasta = out_align

    info = {'ifile':filename,
            'ofile':tfasta,
            'maxmem':MAX_MEM,
            'maxhours':1}
    cmd = 'muscle -in %(ifile)s -out %(ofile)s -maxmb %(maxmem)i -maxhours %(maxhours)i'
    args = shlex.split(cmd % info)
    call(args)


def fasta2aln(in_file, out_file):
    """Converts fasta-alignment files into tab-delimited alignments."""
    with open(out_file, 'w') as handle:
        for name, seq in fasta_iter(in_file):
            handle.write('%s\t%s\n' % (name, seq))

def crazy_iter(source_lim, target_lim, widths, last_items = None):
    """An iterator which returns ranges to check for alignments.

    This iterator has a strange return format which tries to take advantage
    of memcached items. It also takes care of restarting incomplete results.

    Arguements:
    source_lim -- The length of the source sequence
    target_lim -- The length of the target sequence
    widths -- The desired widths to check as a LIST ie. [1,2,3,4,5]

    Kwargs:
    last_items -- a tuple of (source-width, target-width, 
                  source-start, target-start) which indicates where
                  to resume iteration.
    
    Yields:
    (source-width, target-width, source-start, target-start) tuples

"""
    
    def pred(args):
        return not all([x == y for x,y in zip(args, last_items)])

    if last_items is not None:
        swidths = [x for x in widths if x >= last_items[0]]
    else:
        swidths = widths
    
    it = product(swidths, widths, range(*source_lim), range(*target_lim))
    if last_items is not None:
        it = dropwhile(pred, it)

    for sw, tw, ss, ts in it:
        if not(sw+ss > source_lim[-1] or tw+ts > target_lim[-1]):
            yield sw, tw, ss, ts

def getOverlap(a, b):
    """Calculates the amount of overlap between two ranges."""
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    
def get_last(iterable):
    """Finds the last completed tuple.
    
    Used for restarting an un-finished analysis. It iterates through the 
    results file and finds the last completed analysis tuple.
    
    Arguements:
    iterable -- A file-like iterable that was created by Prediction Analysis

    Returns:
    Last completed tuple: 
    (source-width, target-width, source-start, target-start)"""

    if tuple(iterable.fieldnames) != tuple(LINK_FIELDS):
        return None

    def process(iterable):
        while True:        
            try:
                row = iterable.next()
            except StopIteration:
                break
            except:
                continue
            try:
                ss = int(row['Source-Start'])
                ts = int(row['Target-Start'])
                sw = int(row['Source-End'])-ss
                tw = int(row['Target-End'])-ts
                yield (sw, tw, ss, ts)
            except:
                pass

    l = deque(process(iterable), maxlen = 1)
    try:
        row = l.pop()
    except IndexError:
        return None

    return row
    

def check_headers(filename):
    try:
        with open(filename) as handle:
            headers = handle.next().strip().split('\t')
            return headers == LINK_FIELDS
    except StopIteration:
        return False

def get_corrected_mutual_info(signal1, signal2, **kwargs):
    return kwargs['Mutual-Info'] - kwargs['Null-Mutual-Info']


def PredictionAnalysis(align1, align2, outfile, cons_cut = 0.99, **kwargs):

    if not os.path.exists(outfile):
        mode = 'w'
        last = None
    elif check_headers(outfile):
        mode = 'w'
        last = None
    else:
        mode = 'a'
        with open(outfile) as handle:
            iterable = csv.DictReader(handle, delimiter = '\t')
            last = get_last(iterable)

    a1 = Alignment.alignment_from_file(align1)
    a2 = Alignment.alignment_from_file(align2)

    sprot = prot_from_path(align1)
    tprot = prot_from_path(align2)

    defaults = dict(zip(LINK_FIELDS, [None]*len(LINK_FIELDS)))
    defaults['Source-Prot']=sprot
    defaults['Target-Prot']=tprot
    defaults.pop('Source-Start')
    defaults.pop('Source-End')
    defaults.pop('Target-Start')
    defaults.pop('Target-End')

    calculator = LinkCalculator()
    rmheaders = dict((head, None) for head in calculator.get_fields())

    headers = sorted(set(a1.seqs.keys()) & set(a2.seqs.keys()))

    iterable = product(izip(count(), a1.iterate_columns(headers)),
                        izip(count(), a2.iterate_columns(headers)))
    if last:
        fiterable = dropwhile(iterable, lambda x:x[0][0]<last[2] and x[1][0]<last[3])
    else:
        fiterable = iterable

    ohandle = open(outfile, mode)
    writer = csv.DictWriter(ohandle, LINK_FIELDS, delimiter = '\t')
    if mode == 'w':
        writer.writerow(dict(zip(LINK_FIELDS, LINK_FIELDS)))

    for (ind1, seq1), (ind2, seq2) in fiterable:

        cseq1 = ''
        cseq2 = ''
        for s1, s2 in zip(seq1, seq2):
            if s1 != '-' and s2 != '-':
                cseq1 += s1
                cseq2 += s2

        if not cseq1 or not cseq2:
            continue

        c1 = make_counts(cseq1)
        c2 = make_counts(cseq2)

        row = dict()
        row.update(defaults)
        row['Source-Start'] = ind1
        row['Source-End'] = ind1+1
        row['Target-Start'] = ind2
        row['Target-End'] = ind2+1

        row['Source-Cons'] = max(x/len(cseq1) for x in c1.values())
        row['Target-Cons'] = max(x/len(cseq2) for x in c2.values())
        if row['Source-Cons'] == 1 or row['Target-Cons'] == 1:
            writer.writerow(row)
            continue

        row['Source-Entropy'] = calculate_entropy(cseq1)
        row['Target-Entropy'] = calculate_entropy(cseq2)

        for field, val in calculator.calculate_all(cseq1, cseq2):
            row[field] = val

        mappings = prediction_mapping(cseq1, cseq2)
        score = sum([z for _, _, z in mappings])/len(cseq1)
        row['Total-Score'] = score
        source, target, val = mappings[0]
        row.update({'Source-Seq':source,
                    'Target-Seq':target,
                    'Correct-Num':val,
                    'Total-Num':c1[source],
                    'This-Score': val/c1[source]})
        writer.writerow(row)
        row.update(rmheaders)

        for source, target, val in mappings[1:]:
            row.update({'Source-Seq':source,
                        'Target-Seq':target,
                        'Correct-Num':val,
                        'Total-Num':c1[source],
                        'This-Score': val/c1[source]})
            writer.writerow(row)
        ohandle.flush()
        os.fsync(ohandle.fileno())


def OldPredictionAnalysis(align1, align2, outfile, widths = range(1,5), same = False, mode = 'a', cons_cut = 0.98, calc_pval = False, short_linkage_format = False):
    """Analyzes the linkages between 2 alignments.
    
    A controller function which calculates the Linkage between columns in 
    the alignment files. This function takes care of opening/loading 
    alignments, iterating through the columns, calculating linkages, 
    and writing to an output file.

    Arguements:
    align1 -- Path to source alignment file.
    align2 -- Path to target alignment file.
    outfile -- Path to the results file.
    
    Kwargs:
    widths -- The desired column widths to check. Default: range(1,5)
    same -- A boolean indicating whether these files are the same protein.
            Default: False
    mode -- Which mode to open the results file. Default: 'a'
    cons_cut -- The conservation cutoff to use for ignoring columns. 
                Default: 0.8
    calc_pval -- A boolean indicating whether to calculate p-values for 
                 linkages. Default: False
    
    Returns:
    None"""
    
    def get_signals(align1, align2, widths, same, last):
        for sw, tw, ss, ts in crazy_iter([0, align1.width], [0, align2.width], widths, last_items = last):
            if same and getOverlap((ss, ss+sw), (ts, ts+tw)) > 0:
                continue
            if (ss, ss+sw) not in source_skip and (ts, ts+tw) not in target_skip:
                a1 = align1.get_slice(ss, ss+sw)
                a2 = align2.get_slice(ts, ts+tw)
                if a1 is None:
                    source_skip.add((ss, ss+sw))
                if a2 is None:
                    target_skip.add((ts, ts+tw))
                if a1 is None or a2 is None:
                    continue
                    
                over = set(a1.seqs.keys()) & set(a2.seqs.keys())
                if len(over) > 5:
                    yield a1, a2, sorted(over), {'Source-Start':ss, 'Source-End':ss+sw,
                                                    'Target-Start':ts, 'Target-End':ts+tw}
                else:
                    yield None, None, None, {'Source-Start':ss, 'Source-End':ss+sw,
                                                    'Target-Start':ts, 'Target-End':ts+tw}
                if len(a1.seqs) < 10:
                    source_skip.add((ss, ss+sw))
                if len(a2.seqs) < 10:
                    target_skip.add((ts, ts+tw))
                    
                    

    print 'widths!', widths

    line_count = 0
    a1 = Alignment.alignment_from_file(align1)
    a2 = Alignment.alignment_from_file(align2)

    sprot = prot_from_path(align1)
    tprot = prot_from_path(align2)

    defaults = dict(zip(LINK_FIELDS, [None]*len(LINK_FIELDS)))
    defaults['Source-Prot']=sprot
    defaults['Target-Prot']=tprot
    defaults.pop('Source-Start')
    defaults.pop('Source-End')
    defaults.pop('Target-Start')
    defaults.pop('Target-End')

    source_skip = set()
    target_skip = set()

    calculator = LinkCalculator()
    
    if mode == 'a' and os.path.exists(outfile):
        print 'trying to get last line!'
        with open(outfile) as handle:
            last = get_last(csv.DictReader(handle, delimiter = '\t'))
    else:
        last = None

    with open(outfile, mode) as handle:
        handle.write('\t'.join(LINK_FIELDS)+'\n')
        writer = csv.DictWriter(handle, fieldnames = LINK_FIELDS, 
                                delimiter = '\t')
        for slice1, slice2, seqs, loc in get_signals(a1, a2, widths, same, last):
            loc.update(defaults)            
            if slice1 is None:
                loc.update({'Source-Seq':None,
                            'Target-Seq':None,
                            'Correct-Num':'too few',
                            'Total-Num':'too few',
                            'This-Score': 0})
                #print 'few %(Source-Start)i, %(Source-End)i, %(Target-Start)i, %(Target-Start)i' % loc
                if not short_linkage_format:                
                    writer.writerow(loc)
                continue
            if not loc['Target-Start'] % 10:
                print '%(Source-Prot)s,%(Target-Prot)s,%(Source-Start)i,%(Target-Start)i' % loc
            s1, m1 = slice1.get_signal(seqs)
            s2, m2 = slice2.get_signal(seqs)

            #create reverse mappings
            rm1 = dict([(y,x) for x,y in m1.items()])
            rm2 = dict([(y,x) for x,y in m2.items()])

            seq1 = ''.join(rm1[s].upper() for s in s1)
            seq2 = ''.join(rm2[s].upper() for s in s2)

            #create count dictionary
            c1 = make_counts(s1)
            c2 = make_counts(s2)
            
            loc['Source-Cons'] = max(x/len(s1) for x in c1.values())
            loc['Target-Cons'] = max(x/len(s2) for x in c2.values())

            if loc['Source-Cons'] > cons_cut or loc['Target-Cons'] > cons_cut:
                loc.update({'Correct-Num':'too conserved',
                            'Total-Num':'too conserved'})
                if not short_linkage_format:
                    writer.writerow(loc)
                if loc['Source-Cons'] > cons_cut:
                    source_skip.add((loc['Source-Start'], loc['Source-End']))
                if loc['Target-Cons'] > cons_cut:
                    target_skip.add((loc['Target-Start'], loc['Target-End']))
                continue


            mappings = prediction_mapping(tuple(s1), tuple(s2))
            score = sum([z for _, _, z in mappings])/len(s1)
            loc['Total-Score'] = score
            for field, val in calculator.calculate_all(seq1, seq2):
                loc[field] = val

            #print '%(Source-Start)i, %(Source-End)i, %(Target-Start)i, %(Target-Start)i, %(Total-Score)f' % loc
            line_count += 1
            for source, target, val in mappings:
                loc.update({'Source-Seq':rm1[source],
                                'Target-Seq':rm2[target],
                                'Correct-Num':val,
                                'Total-Num':c1[source],
                                'This-Score': val/c1[source]})
                writer.writerow(loc)
            handle.flush()
            os.fsync(handle.fileno())
    return line_count











