import sqlite3
from types import ListType, TupleType
from itertools import product
import os.path
import hashlib



def hexhash(path):
    """Returns the md5 hexhash of the file at the path."""

    m = hashlib.md5()
    with open(path) as handle:
        for line in handle:
            m.update(line)
    return m.hexdigest()


def need_to_do(ifiles, ofiles, *args, **kwargs):
    """Checks whether the files need to be updated."""

    ohashes = {}
    if type(ofiles) == ListType or type(ofiles) == TupleType:
        for f in ofiles:
            if not os.path.exists(f):
                return True, 'Missing file: %s' % f
            else:
                ohashes[f] = hexhash(f)
    else:
        if not os.path.exists(ofiles):
            return True, 'Missing file: %s' % f
        else:
            ohashes[ofiles] = hexhash(ofiles)

    if type(ifiles) == ListType or type(ifiles) == TupleType:
        ihashes = dict([(x, hexhash(x)) for x in ifiles])
    else:
        ihashes = dict()
        ihashes[ifiles] = hexhash(ifiles)

    con = sqlite3.connect('filedata.sql')
    stri = "select * from dep where spath=? and shash=? and dpath=? and dhash=?"
    iterable = product(ihashes.iteritems(), ohashes.iteritems())
    for (spath, shash), (dpath, dhash) in iterable:
        rows = con.execute(stri, (spath, shash, dpath, dhash))
        if rows.fetchone() is None:
            return True, 'Files out of date!'
    
    return False, 'All files up to date!'
