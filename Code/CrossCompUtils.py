import csv
from functools import partial
from operator import itemgetter
from itertools import *

try:
    from collections import Counter
except ImportError:
    from collections import defaultdict
    Counter = partial(defaultdict, int)
        


def linkage_iter(link_file, organism):
    """Iterates single linkages for a file"""
    
    with open(link_file) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            row['Organism'] = organism
            yield row


def compare_linkages(link_files, org_names, outfiles):
    """Compares the linkages between many files and outputs a summary."""
    
    def get_prot_overlaps(input_dict):
        """Yields any rows of the corresponding protiens or returns None"""
        
        keyfun = itemgetter('Source-Prot', 'Target-Prot')
        
        grouped = {}
        for key, item in input_dict.items():
            grouped[key] = groupby(item, keyfun)
    
        current_items = {}
        for key, val in grouped.items():
            current_items[key] = val.next()
    
        while current_items:
            this_item = min((item for item, _ in current_items.itervalues()))
            outdict = {}
            print this_item
            for key, (item, rows) in current_items.items():
                if item == this_item:
                    outdict[key] = list(rows)
                    try:
                        current_items[key] = grouped[key].next()
                    except StopIteration:
                        current_items.pop(key)
                        grouped.pop(key)
            yield this_item, outdict
    
    
    simple_fields = ['Source-Prot', 'Target-Prot'] + sorted(org_names)
    complex_fields = ('Organism', 'Source-Prot', 'Target-Prot','Source-Start',
                'Source-End','Target-Start','Target-End', 'Total-Score')
    
    simple_writer = csv.DictWriter(open(outfiles[0], 'w'), simple_fields, 
                                    delimiter = '\t', extrasaction = 'ignore')
    complex_writer = csv.DictWriter(open(outfiles[1], 'w'), complex_fields, 
                                    delimiter = '\t', extrasaction = 'ignore')
    simple_writer.writerow(dict(zip(simple_fields, simple_fields)))
    complex_writer.writerow(dict(zip(complex_fields, complex_fields)))
    
    iterable_dict = {}
    for org, fname in zip(org_names, link_files):
        iterable_dict[org] = linkage_iter(fname, org)
    
    for prots, outdict in get_prot_overlaps(iterable_dict):
        counter = Counter()
        print prots
        for org, rows in sorted(outdict.items()):
            for row in rows:
                complex_writer.writerow(row)
                counter[org] += 1
        counter['Source-Prot'], counter['Target-Prot'] = prots
        simple_writer.writerow(counter)
        
        
        
        
        
        
        
        
