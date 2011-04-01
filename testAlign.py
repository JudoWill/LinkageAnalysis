import nose
from Code import AlignUtils
import os, os.path

ALNDATA = (('key1', 'AACAA'), ('key2', 'A-CAA'), ('key3', 'TA-AA'))
ALNFILE = os.path.join(os.environ.get('TMPDIR', '/tmp/'), 'tmp.aln')

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

@nose.tools.with_setup(aln_write_fun, aln_delete)        
def test_getting_none_alignment_slice():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    new_aln = aln.get_slice(3,5)
    nose.tools.eq_(new_aln, None)
        
@nose.tools.with_setup(aln_write_fun, aln_delete)        
def test_signal_function():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    wanted_keys = ('key1', 'key2', 'key3')
    signal, mapping = aln.get_signal(wanted_keys)
    
    nose.tools.eq_(signal, range(1,len(wanted_keys)+1))
    rmapping = {}
    for (_, seq), num in zip(ALNDATA, range(1,len(wanted_keys)+1)):
        rmapping[seq] = num
    nose.tools.eq_(mapping, rmapping)

@nose.tools.with_setup(aln_write_fun, aln_delete)    
def test_write_fasta():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    tfile = os.path.join(os.environ.get('TMPDIR', '/tmp/'), 'ntmp.aln')
    try:
        aln.write_fasta(tfile)
        with open(tfile) as handle:
            indata = handle.read()
        
        for key, seq in ALNDATA:
            assert '>%s\n' % key in indata
            assert seq+'\n' in indata
        
    finally:
        os.remove(tfile)
        
@nose.tools.with_setup(aln_write_fun, aln_delete)    
def test_write_phylip():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)        
    aln.seqs['reallylongkeynamethatshouldbeshortened'] = 'ATCTG'

    tfile = os.path.join(os.environ.get('TMPDIR', '/tmp/'), 'ntmp.aln')
    try:
        aln.write_phylip(tfile)
        with open(tfile) as handle:
            fline = handle.next()
            indata = ''.join([line for line in handle])
        cline = '%i %i\n' % (len(aln.seqs), aln.width)
        nose.tools.eq_(fline, cline)
        
        for key, seq in aln.seqs.items():
            assert '%s' % key[:10] in indata
            assert seq+'\n' in indata
        
    finally:
        os.remove(tfile)
        
@nose.tools.with_setup(aln_write_fun, aln_delete)    
def test_get_consensus():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    seq = aln.get_consensus()
    nose.tools.eq_(seq, 'AACAA')       

@nose.tools.with_setup(aln_write_fun, aln_delete)    
def test_get_reference_numbering():
    
    aln = AlignUtils.Alignment.alignment_from_file(ALNFILE)
    
    (dest_nums, align_nums) = aln.convert_numbering('key2')
    nose.tools.eq_(dest_nums, [0,0,1,2,3])
    nose.tools.eq_(align_nums, [0,2,3,4])
