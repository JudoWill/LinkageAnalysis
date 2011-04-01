import nose
from Code import AlignUtils
import os, os.path

ALNDATA = (('key1', 'AACAA'), ('key2', 'A-CAA'), ('key3', 'TA-AA'))
ALNFILE = os.path.join(os.environ['TMPDIR'], 'tmp.aln')

def aln_write_fun():
    
    with open(ALNFILE, 'w') as handle:
        for key, seq in ALNDATA:
            handle.write('%s\t%s\n' % (key,seq))
    
def aln_delete():
    
    try:
        os.remove(ALNFILE)
    except IOError:
        pass

def test_load_alignment():
    
    aln = AlignUtils.Alignment()

@nose.tools.with_setup(aln_write_fun, aln_delete)
def test_load_alignment_from_file():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    nose.tools.eq_(aln.width, 5)
    
    for key, seq in ALNDATA:
        nose.tools.eq_(aln.seqs[key], seq)

@nose.tools.with_setup(aln_write_fun, aln_delete)            
def test_join_alignments_joining_sequences():
    
    aln1 = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    aln2 = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    aln1.append_alignment(aln2)
    nose.tools.eq_(aln1.width, 10)
    nose.tools.eq_(len(aln1.seqs), len(ALNDATA))
    
    for key, seq in ALNDATA:
        nose.tools.eq_(aln1.seqs[key], seq+seq)    
        
@nose.tools.with_setup(aln_write_fun, aln_delete)            
def test_join_alignments_limiting_sequences():
    
    aln1 = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    aln2 = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    aln2.seqs['junkkey'] = 'ATCAA'
    
    
    aln1.append_alignment(aln2)
    nose.tools.eq_(aln1.width, 10)
    
    for key, seq in ALNDATA:
        nose.tools.eq_(aln1.seqs[key], seq+seq)
        
    nose.tools.eq_(len(aln1.seqs), len(ALNDATA))
        
@nose.tools.with_setup(aln_write_fun, aln_delete)        
def test_getting_alignment_slice():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    new_aln = aln.get_slice(1,3)
    print new_aln.seqs
    for key, seq in ALNDATA:
        nose.tools.eq_(new_aln.seqs[key], seq[1:3])      
        
        
        
        
        