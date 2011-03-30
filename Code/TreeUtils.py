from AlignUtils import Alignment
from GeneralUtils import pushd
from subprocess import call
import shlex, os.path, os
from StringIO import StringIO
from dendropy import Tree




def run_phylip(direc, progtype, input_args = ['y'], capture_output = False,
                clean_direc = True):
    """Runs the phylip suite of programs"""

    with pushd(direc):
        with open('input', 'w') as handle:
            handle.write('\n'.join(input_args+['']))
        inbuf = open('input')
        if progtype == 'proml':
            check_files = ('infile', input_args[0])
            rmfiles = ('outfile', 'outtree')
            cmd = 'phylip proml'
        elif progtype == 'consense':
            check_files = ('intree', input_args[0])
            rmfiles = ('outfile', 'outtree')
            cmd = 'phylip consense'
        else:
            raise TypeError, 'Unrecognized phylip program: %s' % progtype

        if clean_direc:
            for f in rmfiles:
                try:
                    os.remove(f)
                except OSError:
                    pass

        if not any([os.path.exists(x) for x in check_files]):
            raise OSError, 'Could not find any input files!'
        
        if capture_output:
            out = open('screenout', 'w')
        else:
            out = open(os.devnull, 'w')

        args = shlex.split(cmd)
        call(args, stdin = inbuf, stdout = out)

def group_sequences(tree_file, cutoff):
    """Loads a tree to find sequences which can be merged."""

    with open(tree_file) as handle:
        tree = Tree.get_from_string(handle.read(), schema = 'newick')

    groups = []
    done = set()
    for t in tree.level_order_node_iter():
        if t.edge.length >= cutoff:
            lens = [x.edge.length >= cutoff for x in t.level_order_iter() if not x.is_leaf()]
            if all(lens) or not any(lens) or t.is_leaf():
                leafs = set([str(x.taxon) for x in t.leaf_nodes()])
                groups.append(leafs-done)
                done |= leafs

    groups = [x for x in groups if x]
    return groups


def merge_sequences(ifile, ofile, tree_file, cutoff = 70):
    """Merges sequence alignments based on a tree"""

    def make_consensus(tups):
        aln = Alignment()
        aln.seqs = dict(tups)
        aln.width = len(tups[0][1])
        return aln.get_consensus()

    groups = [x for x in group_sequences(tree_file, cutoff) if len(x) > 1]

    original = Alignment.alignment_from_file(ifile)        
    for group in groups:
        tup = []
        for key in group:
            tup.append((key, original.seqs.pop(key)))
        nkey = '-'.join(sorted(group))
        original.seqs[nkey] = make_consensus(tup)
    original.write_aln(ofile)
    original.write_fasta(ofile+'.fasta')
            























