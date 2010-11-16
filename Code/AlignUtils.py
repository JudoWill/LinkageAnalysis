from __future__ import division
from collections import deque, defaultdict
from GeneralUtils import *
from subprocess import call
import shlex
from itertools import groupby
from math import log
from memorised.decorators import memorise
from random import shuffle
from operator import itemgetter


class Alignment():
    
    def __init__(self):
        self.seqs = {}
    
    @staticmethod
    def alignment_from_file(filename):
        align = Alignment()        
        with open(filename) as handle:
            for line in handle:
                parts = line.strip().split('\t')
                align.seqs[parts[0]] = align.seqs[parts[1]]
        return align

    def get_slice(self, start, stop):
        
        align = Alignment()
        for name, seq in self.seqs:
            align.seqs[name] = seq[start:stop]
        return align

    def get_signal(self, seq_names):
        
        signal = []
        count = 0
        seq_dict = {}
        for name in seq_names:
            seq = self.seqs[name]
            signal.append(seq_dict.setdefault(seq, len(seq_dict)+1))

        return signal

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
        if calculate_mutual_info(signal1, shuffle(signal2)) > rmut:
            num_greater += 1

    return num_greater / num_reps

@memorise()
def prediction_mapping(signal1, signal2):
    print signal1
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
    
    join_fasta(filenames, out_fasta)
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

def join_alignments(out_align, *aligns):
    
    base_align = aligns[0]
    for align in aligns[1:]:
        info = {
                'align1':base_align,
                'align2':align,
                'out':out_align
                }
        cmd = 'clustalw -PROFILE1=%(align1)s -PROFILE2=%(align2)s -QUICKTREE -OUTFILE=%(out)s -ENDGAPS -PROFILE -QUIET'
        args = shlex.split(cmd % info)
        call(args)
        base_align = out_align

def convert_alignment(clustal_v, modified_v):
    
    grouper = lambda x: x.startswith(' ') or x.startswith('\n')
    seqs = defaultdict(str)    
    
    with open(clustal_v) as clustalHandle:
        clustalHandle.next()
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


