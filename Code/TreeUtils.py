from AlignUtils import Alignment
from GeneralUtils import pushd
from subprocess import call
import shlex, os.path



def run_proml(direc, input_args = ['Y'], capture_output = False, 
                clean_direc = True):
    """Runs the phylip proml function."""

   with pushd(direc):
        #make sure file is there
        fpresent = os.path.exists('infile') or os.path.exists(input_args[0].strip())
        if not fpresent:
            rasie IOError, 'Could not find input file!'

        with open('input', 'w') as handle:
            handle.write('\n'.join(input_args))
        
        if clean_direc:
            outs = ('outfile', 'outtree')
            for f in outs:
                try:
                    os.remove(f)
                except OSError:
                    pass

        cmd = 'phylip proml < input'
        if capture_output:
            cmd += ' > screenout'
        else:
            cmd += ' > /dev/null'

        args = shlex.split(cmd)
        call(args)
