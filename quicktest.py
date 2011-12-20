import os, os.path, csv
from itertools import *
from operator import itemgetter
from Code.AlignUtils import Alignment
from Code.StatUtils import linkage_pval, check_linkage_file

if __name__ == '__main__':

    outfile = 'StaphData/SigResults.txt'
    link_dir = 'StaphData/LinkageResults/'
    fields = ('Source-Prot', 'Target-Prot', 'Source-Start','Source-End',
                'Target-Start','Target-End', 'Total-Num',
                'This-Score', 'Total-Score', 'p-val')
    
    files = sorted([x for x in os.listdir(link_dir) if x.endswith('.res')])
    if os.path.exists(outfile):
        t = itemgetter(*fields[:3])
        arows = []
        with open(outfile) as handle:
            n = []
            for key, rows in groupby(csv.DictReader(handle, delimiter = '\t'), t):
                arows += n
                n = list(rows)
            lspot = n[0]['Source-Prot'] + '--' + n[0]['Target-Prot'] + '.res'
        files = list(ifilter(lambda x: x > lspot, files))
        handle = open(outfile, 'w')
        writer = csv.DictWriter(handle, fields, delimiter = '\t', extrasaction = 'ignore')
        writer.writerow(dict(zip(fields, fields)))
        writer.writerows(arows)
        
    else:
        handle = open(outfile, 'w')
        writer = csv.DictWriter(handle, fields, delimiter = '\t', extrasaction = 'ignore')
        writer.writerow(dict(zip(fields, fields)))
        
        
    tstr = '%(Source-Prot)s\t%(Target-Prot)s\t%(Source-Start)s\t%(Target-Start)s\t%(p-val)f'
    for f in files:
        print 'testing', f
        for row in check_linkage_file(os.path.join(link_dir, f), 
                                        'StaphData/LANLSequences/Alignments/'):
            if row['p-val'] < 0.1:
                
                print tstr % row
            writer.writerow(row)
    
    