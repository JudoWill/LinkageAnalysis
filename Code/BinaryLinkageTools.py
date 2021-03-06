from struct import Struct
from LinkFields import LINK_FIELDS
from operator import itemgetter, methodcaller, attrgetter
from collections import namedtuple
from itertools import chain, tee, islice
from functools import partial
import csv
from fileinput import FileInput


class BinaryLinkageResults(list):

    
    def __init__(self, fname, strsize = 8):

        self.f = open(fname, 'rb')
        self.fname = fname
        self.strsize = strsize
        self.struct, self.mapping_list = BinaryLinkageResults.generate_struct(LINK_FIELDS,8)
        self.f.seek(0, 2)
        self.len = self.f.tell()/self.struct.size
        self.bsize = self.struct.size
        self.key = itemgetter(*range(6))
        self.fields = [x.replace('-', '_') for x in LINK_FIELDS]
        self.tuple_class = namedtuple('record', self.fields)

    
    @classmethod
    def generate_struct(fields, strsize):
    
        return_next = methodcaller('next')
        mapping = (('Start','i', return_next), ('End','i', return_next), ('Num', 'i', return_next),
            ('Prot', '%ic'%strsize, lambda x: ''.join(islice(x, strsize)).strip()),
            ('Score', 'f', return_next), ('Cons', 'f', return_next), 
            ('Info', 'f', return_next), ('Dist', 'f', return_next),
            ('Seq', 'c', return_next))

        mapping_list = []
        fmtstr = ''
        for field in LINK_FIELDS:
            for key, fmt, func in mapping:
                if key in field:
                    fmtstr += fmt
                    mapping_list.append(func)
                    break

        return Struct(fmtstr), mapping_list


    @classmethod
    def convert_files_to_binary(infiles, outfile):
        
        infiles = sorted(infiles)
        
        struct, mapping_list = BinaryLinkageResults.generate_struct(LINK_FIELDS, 8)
        
        handle = FileInput(infiles)
        reader = csv.reader(handle, delimiter = '\t')
        grouper = itemgetter(*range(6)) #group by source-prot through target-end
        buf = open(outfile, 'wb')
        count = 0
        for ind, (key, rows) in enumerate(groupby(reader, key = grouper)):
            trow = rows.next()
            converted = [func(field) for func, field in zip(mapping_list, trow)]
            if all(x is not None for x in converted):
                arglist = list(chain(converted[0], converted[1], converted[2:]))
                data = StructClass.pack(*arglist)
                buf.write(data)
                count += 1
                if count % 10000 == 0:
                    print count, key



    def __len__(self):
        return self.len

    def iterate_records(self):
        
        tmphandle = open(self.fname, 'rb') #so its not interupted by other things
        sz = self.bsize #for faster lookups
        struct = self.struct #for faster lookups
        convfun = self._convert_item #for faster lookups
        for recnum in xrange(self.len):
            buf = tmphandle.read(sz)
            vals = struct.unpack_from(buf)
            yield convfun(vals)
    
    def _convert_item(self, unpacked_data):
        iterable = iter(unpacked_data)
        
        return self.tuple_class(*[func(iterable) for func in self.mapping_list])

    
    def __getitem__(self,i):
        if i >= 0 or i < self.len:
            self.f.seek(i*self.bsize)
            buf = self.f.read(self.bsize)
            vals = self.struct.unpack_from(buf)
            item = self._convert_item(vals)
            return item
        else:
            raise ValueError

    def _generate_reverse_item(self,item):
        """Generates a reverse item for use as a lookup for pairwise iteractions."""
        
        mapping = dict(zip(self.fields, [None]*len(self.fields)))      
        mapping.update({'Source_Prot':item.Target_Prot, 
                        'Source_Start':item.Target_Start, 
                        'Source_End':item.Target_End, 
                        'Target_Prot':item.Source_Prot, 
                        'Target_Start':item.Source_Start, 
                        'Target_End':item.Source_End})
        
        return self.tuple_class(**mapping)
        
    def get_reverse_item(self, item, ind = 0):
        """Returns the pairwise item for the one provided.

        Raises ValueError if it cannot find it.
        """
        
        rev_item = self._generate_reverse_item(item)
        rev_index = self.index(rev_item, lo = ind)
        return self[rev_index]
        
    def get_pairwise_records(self):

        seen = set()
        
        for ind, record in enumerate(self.iterate_records()):
            if any(x == -1 for x in record) or record.Source_Start == record.Target_Start:
                continue
            try:
                rev_item = self.get_reverse_item(record, ind = ind)
            except ValueError:
                continue
            if rev_item not in seen: #we've reached the middle
                yield record, rev_item
                seen.add(rev_item)


    
    def index(self, x, key = None, lo = 0, hi = None):
        """Locate the leftmost value exactly equal to x"""      
        if key is None:
            key = self.key
        i = bisect_left(self, x, key = key, lo = lo, hi = hi)
        
        if i != len(self) and key(self[i]) == key(x):
            return i
        raise ValueError
        

    




def bisect_left(a, x, lo=0, hi=None, key = lambda x:x):
    """Return the index where to insert item x in list a, assuming a is sorted.

    The return value i is such that all e in a[:i] have e < x, and all e in
    a[i:] have e >= x.  So if x already appears in the list, a.insert(x) will
    insert just before the leftmost x already there.

    Optional args lo (default 0) and hi (default len(a)) bound the
    slice of a to be searched.
    """
	
    if lo < 0:
        raise ValueError('lo must be non-negative')
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        if key(a[mid]) < key(x): lo = mid+1
        else: hi = mid
    return lo



def safefloat(string):
    try:
        return float(string)
    except ValueError:
        if len(string) == 0:
            return 0.0
        else:
            return None

def safeint(string):
    try:
        return int(string.strip())
    except ValueError:
        if len(string) == 0:
            return 0
        else:
            return None
