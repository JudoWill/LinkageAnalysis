from __future__ import division
from collections import deque, defaultdict
from GeneralUtils import *
from subprocess import call
import shlex
from itertools import groupby, product, dropwhile, imap
from math import log
from random import shuffle, sample
from operator import itemgetter, ne, eq
import csv
from functools import partial

try:
    from memorised.decorators import memorise
except ImportError:
    class memorise(object):
        def __init__(self, func):
            self.func = func
        def __call__(self, *args, **kwargs):
            return self.func(*args, **kwargs)
        

class Alignment():
    
    def __init__(self):
        self.seqs = {}
        self.width = None
    
    @staticmethod
    def alignment_from_file(filename):
        align = Alignment()        
        with open(filename) as handle:
            for line in handle:
                parts = line.strip().split('\t')
                align.seqs[parts[0]] = parts[1]
                align.width = len(parts[1])
        return align

    def get_slice(self, start, stop):
        
        ngapfun = partial(ne, '-')
        gapfun = partial(eq, '-')
        align = Alignment()
        for name, seq in self.seqs.items():
            if any(imap(ngapfun, seq[:stop])) and not all(imap(gapfun, seq[start:stop])):
                align.seqs[name] = seq[start:stop]
        return align

    def get_signal(self, seq_names):
        
        signal = []
        count = 0
        seq_dict = {}
        for name in seq_names:
            seq = self.seqs[name]
            signal.append(seq_dict.setdefault(seq, len(seq_dict)+1))

        return signal, seq_dict

@memorise()
def calculate_mutual_info(signal1, signal2):
    
    def count2prob(d, num):
        for key, val in d.items():
            d[key] = val/num
        return d


    overlap = defaultdict(int)
    num_items = len(signal1)
    signal1_hist = defaultdict(int)
    signal2_hist = defaultdict(int)

    for s1, s2 in zip(signal1, signal2):
        overlap[(s1, s2)] += 1
        signal1_hist[s1] += 1
        signal2_hist[s2] += 1

    mut_info = float()
    overlap_prob = count2prob(overlap, num_items)
    signal1_prob = count2prob(signal1_hist, num_items)
    signal2_prob = count2prob(signal2_hist, num_items)

    for (s1, s2), count in overlap.items():
        mut_info += overlap_prob[(s1, s2)]*log(overlap_prob[(s1, s2)]/(signal1_prob[s1]*signal2_prob[s2]))
        
    return mut_info

@memorise()
def get_mutual_info_pval(signal1, signal2, num_reps = 5000):
    
    rmut = calculate_mutual_info(signal1, signal2)

    num_greater = 0
    
    for i in xrange(num_reps):
        
        r = calculate_mutual_info(signal1, sample(signal2, len(signal2)))
        print r, rmut
        if r > rmut:
            num_greater += 1

    return num_greater / num_reps

@memorise()
def prediction_mapping(signal1, signal2):
    print 'calculating!'
    counts = defaultdict(int)
    for s1, s2 in zip(signal1, signal2):
        counts[(s1, s2)] += 1

    (s1, s2), val = max(counts.items(), key = itemgetter(1))
    if val == 1:
        return [(x, y, 1) for (x,y) in counts.keys()]
    else:
        mapping = (s1, s2, val)
        sn1 = tuple()
        sn2 = tuple()
        for i1, i2 in zip(signal1, signal2):
            if i1 != s1:
                sn1 += (i1,)
                sn2 += (i2,)
        if sn1:
            return [mapping] + prediction_mapping(sn1, sn2)
        else:
            return [mapping]


def run_clustalw(filenames, out_fasta, out_tree, out_align):
    
    join_fasta(filenames, out_fasta, strip = True)
    info = {   
            'ifile':out_fasta,
            'tree':out_tree,
            'ofile':out_align
            }
    cmd = 'clustalw -INFILE=%(ifile)s -QUICKTREE -NEWTREE=%(tree)s -OUTFILE=%(ofile)s -ENDGAPS -QUIET'
    args = shlex.split(cmd % info)
    call(args)
    
    cmd = 'clustalw -INFILE=%(ifile)s -QUICKTREE -USETREE=%(tree)s -OUTFILE=%(ofile)s -ENDGAPS -QUIET'
    args = shlex.split(cmd % info)
    call(args)

def join_alignments(*aligns):
    
    base_align = aligns[0]
    for align in aligns[1:]:
        info = {
                'align1':base_align,
                'align2':align,
                'out':base_align
                }
        cmd = 'clustalw -PROFILE1=%(align1)s -PROFILE2=%(align2)s -QUICKTREE -OUTFILE=%(out)s -ENDGAPS -PROFILE -QUIET'
        args = shlex.split(cmd % info)
        call(args)

def convert_alignment(clustal_v, modified_v):
    
    grouper = lambda x: x.startswith(' ') or x.startswith('\n')
    seqs = defaultdict(str)    
    
    with open(clustal_v) as clustalHandle:
        try:
            clustalHandle.next()
        except StopIteration:
            return
        for key, rows in groupby(clustalHandle, grouper):
            if not key:
                for row in rows:
                    parts = row.split()
                    seqs[parts[0].strip()] += parts[1].strip()

    num = [len(x) for x in seqs.itervalues()]
    assert all([num[0] == x for x in num])
    
    with open(modified_v, 'w') as handle:
        for name, seq in seqs.iteritems():
            handle.write('%s\t%s\n' % (name, seq))

def crazy_iter(source_lim, target_lim, widths, last_items = None):
    
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



def PredictionAnalysis(align1, align2, outfile, widths = range(1,5)):

    
    def get_signals(align1, align2, widths):
        for sw, tw, ss, ts in crazy_iter([0, align1.width], [0, align2.width], widths):
            if (ss, ss+sw) not in source_skip and (ts, ts+tw) not in target_skip:
                a1 = align1.get_slice(ss, ss+sw)
                a2 = align2.get_slice(ts, ts+tw)
                over = set(a1.seqs.keys()) & set(a2.seqs.keys())
                if len(over) > 300:
                    yield a1, a2, sorted(over), {'Source-Start':ss, 'Source-End':ss+sw,
                                                    'Target-Start':ts, 'Target-End':ts+tw}
                else:
                    yield None, None, None, {'Source-Start':ss, 'Source-End':ss+sw,
                                                    'Target-Start':ts, 'Target-End':ts+tw}
                if len(a1.seqs) < 300:
                    source_skip.add((ss, ss+sw))
                if len(a2.seqs) < 300:
                    target_skip.add((ts, ts+tw))
                    
                    

        
    def make_counts(signal):
        cdict = defaultdict(int)
        for s in signal:
            cdict[s] += 1
        return cdict


    a1 = Alignment.alignment_from_file(align1)
    a2 = Alignment.alignment_from_file(align2)

    fields = ('Source-Start', 'Source-End', 'Target-Start', 'Target-End',
                'Source-Seq', 'Target-Seq', 'Correct-Num', 'Total-Num', 
                'This-Score', 'Total-Score')
    source_skip = set()
    target_skip = set()

    with open(outfile, 'w') as handle:
        handle.write('\t'.join(fields)+'\n')
        writer = csv.DictWriter(handle, fieldnames = fields, 
                                delimiter = '\t')
        for slice1, slice2, seqs, loc in get_signals(a1, a2, widths):
            if slice1 is None:
                loc.update({'Source-Seq':None,
                            'Target-Seq':None,
                            'Correct-Num':'too few',
                            'Total-Num':'too few',
                            'This-Score': 0})
                print 'few %(Source-Start)i, %(Source-End)i, %(Target-Start)i, %(Target-Start)i' % loc
                writer.writerow(loc)
                continue
                          
            s1, m1 = slice1.get_signal(seqs)
            s2, m2 = slice2.get_signal(seqs)
            
            #create reverse mappings
            rm1 = dict([(y,x) for x,y in m1.items()])
            rm2 = dict([(y,x) for x,y in m2.items()])

            #create count dictionary
            c1 = make_counts(s1)
            c2 = make_counts(s2)

            if any([x/len(s1) > 0.8 for x in c1.values()]) or any([x/len(s2) > 0.8 for x in c2.values()]):
                loc.update({'Source-Seq':None,
                                'Target-Seq':None,
                                'Correct-Num':'too conserved',
                                'Total-Num':'too conserved',
                                'This-Score': 0})
                writer.writerow(loc)
                print 'conserved %(Source-Start)i, %(Source-End)i, %(Target-Start)i, %(Target-Start)i' % loc
                if any([x/len(s1) > 0.8 for x in c1.values()]):
                    source_skip.add((loc['Source-Start'], loc['Source-End']))
                if any([x/len(s2) > 0.8 for x in c2.values()]):
                    target_skip.add((loc['Target-Start'], loc['Target-End']))
                print len(s1),c1.values()
                print len(s2),c2.values() 
                continue

            
            mappings = prediction_mapping(tuple(s1), tuple(s2))
            score = sum([z for _, _, z in mappings])/len(s1)
            loc['Total-Score'] = score
            print '%(Source-Start)i, %(Source-End)i, %(Target-Start)i, %(Target-Start)i, %(Total-Score)f' % loc
            for mapping in mappings:
                loc.update({'Source-Seq':rm1[mapping[0]],
                                'Target-Seq':rm2[mapping[1]],
                                'Correct-Num':mapping[2],
                                'Total-Num':c1[mapping[0]],
                                'This-Score': mapping[2]/c1[mapping[0]]})                
                writer.writerow(loc)
            handle.flush()
            os.fsync(handle.fileno())












