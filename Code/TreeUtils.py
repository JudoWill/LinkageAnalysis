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

        


def run_proml(direc, input_args = ['y'], capture_output = False, 
                clean_direc = True):
   with pushd(direc):
        #make sure file is there
        fpresent = os.path.exists('infile') or os.path.exists(input_args[0].strip())
        if not fpresent:
            raise IOError, 'Could not find input file!'

        with open('input', 'w') as handle:
            handle.write('\n'.join(input_args)+'\n')

        if clean_direc:
            outs = ('outfile', 'outtree')
            for f in outs:
                try:
                    os.remove(f)
                except OSError:
                    pass

        cmd = 'phylip proml'
        if capture_output:
            out = open('screenout')
        else:
            out = None

        print cmd
        args = shlex.split(cmd)
        call(args, stdin = open('input'), stdout = out)

def run_consense(direc, input_args = ['Y'], capture_output = False, 
                clean_direc = True):
    """Runs the phylip consense function."""

    with pushd(direc):
        #make sure file is there
        fpresent = os.path.exists('intree') or os.path.exists(input_args[0].strip())
        if not fpresent:
            raise IOError, 'Could not find input file!'

        with open('input', 'w') as handle:
            handle.write('\n'.join(input_args))

        if clean_direc:
            outs = ('outfile', 'outtree')
            for f in outs:
                try:
                    os.remove(f)
                except OSError:
                    pass

        cmd = 'phylip consense < input'
        if capture_output:
            cmd += ' > screenout'
        else:
            cmd += ' > /dev/null'

        args = shlex.split(cmd)
        call(args)

