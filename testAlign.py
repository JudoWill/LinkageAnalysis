import nose
from Code import AlignUtils
import os, os.path


def aln_write_fun():
    
    tfile = os.path.join(os.environ['TMPDIR'], 'tmp.aln')
    indata = (('key1', 'AACAA'), ('key2', 'A-CAA'))
    
    with open(tfile, 'w') as handle:
        for key, seq in indata:
            handle.write('%s\t%s\n' % (key,seq))
    
def aln_delete():
    
    tfile = os.path.join(os.environ['TMPDIR'], 'tmp.aln')
    try:
        os.remove(tfile)
    except IOError:
        pass

def test_load_alignment():
    
    aln = AlignUtils.Alignment()

@nose.tools.with_setup(aln_write_fun, aln_delete)
def test_load_alignment_from_file():
    
    tfile = os.path.join(os.environ['TMPDIR'], 'tmp.aln')
    
    aln = AlignUtils.Alignment.alignment_from_file(tfile)
    nose.tools.eq_(aln.width, 5)
    indata = (('key1', 'AACAA'), ('key2', 'A-CAA'))
    for key, seq in indata:
        nose.tools.eq_(aln.seqs[key], seq)

@nose.tools.with_setup(aln_write_fun, aln_delete)            
def test_join_alignments_joining_sequences():
    
    tfile = os.path.join(os.environ['TMPDIR'], 'tmp.aln')
    aln1 = AlignUtils.Alignment.alignment_from_file(tfile)
    aln2 = AlignUtils.Alignment.alignment_from_file(tfile) 
    indata = (('key1', 'AACAA'), ('key2', 'A-CAA'))
    
    aln1.append_alignment(aln2)
    nose.tools.eq_(aln1.width, 10)
    nose.tools.eq_(len(aln1.seqs), 2)
    
    for key, seq in indata:
        nose.tools.eq_(aln1.seqs[key], seq+seq)    
        
@nose.tools.with_setup(aln_write_fun, aln_delete)            
def test_join_alignments_limiting_sequences():
    
    tfile = os.path.join(os.environ['TMPDIR'], 'tmp.aln')
    aln1 = AlignUtils.Alignment.alignment_from_file(tfile)
    aln2 = AlignUtils.Alignment.alignment_from_file(tfile)
    aln2.seqs['key3'] = 'ATCAA'
    indata = (('key1', 'AACAA'), ('key2', 'A-CAA'))
    
    aln1.append_alignment(aln2)
    nose.tools.eq_(aln1.width, 10)
    
    for key, seq in indata:
        nose.tools.eq_(aln1.seqs[key], seq+seq)
        
    nose.tools.eq_(len(aln1.seqs), 2)
        
        
        
        
        
        
        
        