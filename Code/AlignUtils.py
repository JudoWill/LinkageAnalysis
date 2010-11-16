from collections import deque, defaultdict
from GeneralUtils import *
from subprocess import call
import shlex
from itertools import groupby

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
