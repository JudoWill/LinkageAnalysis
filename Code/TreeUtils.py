from AlignUtils import Alignment
from GeneralUtils import pushd
from subprocess import call
import shlex, os.path, os
from StringIO import StringIO




def run_phylip(direc, progtype, input_args = ['y'], capture_output = False,
                clean_direc = True):
    """Runs the phylip suite of programs"""

    with pushd(direc):
        inbuf = StringIO('\n'.join(input_args+['']))
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
            for f in outs:
                try:
                    os.remove(f)
                except OSError:
                    pass
        
        if capture_output:
            out = open('screenout', 'w')
        else:
            out = open(os.devnull, 'w')

        args = shlex.split(cmd)
        call(args, stdin = inbuf, stdout = out)


