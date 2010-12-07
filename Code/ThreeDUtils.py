from __future__ import division
from copy import deepcopy
from itertools import groupby, product, combinations
from operator import itemgetter
from math import sqrt



conv_dict = {'ala':'A', 'arg':'R', 'asn':'N',
                'asp':'D', 'cys':'C', 'glu':'E',
                'gly':'G', 'his':'H', 'ile':'I',
                'leu':'L', 'lys':'K', 'met':'M',
                'phe':'F', 'pro':'P', 'ser':'S',
                'thr':'T', 'trp':'W', 'tyr':'Y',
                'val':'V', 'gln':'Q'}


def euc_dist(v1, v2):
    return sqrt(sum([(x-y)**2 for x, y in zip(v1, v2)]))


class Structure():
    
    def __init__(self, chain):
        self.chain = chain
        seq = list()
        for num, lines in groupby(chain, itemgetter('resnum')):
            seq.append(conv_dict[lines.next()['residue'].lower()])
        self.seq = ''.join(seq)
        self.pairwise = None

    @staticmethod
    def from_file(filename, chain):
        d = {
            'atom':lambda x: x[2],
            'residue':lambda x: x[3],
            'chain':lambda x: x[4],
            'resnum':lambda x: int(x[5]),
            'X':lambda x: float(x[6]),
            'Y':lambda x: float(x[7]),
            'Z':lambda x: float(x[8]),
            'element':lambda x: x[-1]
                }
        chain_list = list()

        with open(filename) as handle:
            for line in handle:
                if line.startswith('ATOM'):
                    ndict = {}
                    parts = line.split()
                    for key, fun in d.items():
                        ndict[key] = fun(parts)
                    if ndict['chain'] == chain:
                        chain_list.append(deepcopy(ndict))

        return Structure(chain_list)

    def calculate_pairwise(self):
        
        pairwise = {}
        chains = list()
        for resnum, residues in groupby(self.chain, itemgetter('resnum')):
            reslist = list(residues)
            X = sum([x['X'] for x in reslist])/len(reslist)
            Y = sum([x['Y'] for x in reslist])/len(reslist)
            Z = sum([x['Z'] for x in reslist])/len(reslist)
            chains.append((resnum, X, Y, Z))

        for it1, it2 in combinations(chains, 2):
            d = euc_dist(it1[1:], it2[1:])            
            pairwise[(it1[0], it2[0])] = d
            pairwise[(it2[0], it1[0])] = d
        self.pairwise = pairwise

        
