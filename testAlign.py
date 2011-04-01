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
            
        
        
        
        
        
        
        
        
        
        
        
        