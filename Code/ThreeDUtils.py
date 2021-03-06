from __future__ import division
import re
from copy import deepcopy
from itertools import groupby, product, combinations
from operator import itemgetter
from math import sqrt
from nwalign import global_align
from Code.AlignUtils import *
from Code.GeneralUtils import *
from collections import defaultdict
from BeautifulSoup import BeautifulStoneSoup


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
    
    @staticmethod    
    def from_xml(filename, chain):
        
        def get_item(obj, item, conv = str):
            try:
                tmp = obj.find('pdbx:'+item).string
            except AttributeError:
                return None
            try:
                return conv(tmp)
            except TypeError:
                return None
        
        soup = BeautifulStoneSoup(open(filename))
        chain_list = list()
        
        edges = (('pdbx:length_a', 'X'),
                    ('pdbx:length_b', 'Y'),
                    ('pdbx:length_c', 'Z'))
                    
        scales = {}
        for name, edge in edges:
            scales[edge] = float(soup.find(name).contents[0])
        
        for atom in soup.findAll('pdbx:atom_site'):
            if get_item(atom, 'group_pdb') == 'ATOM' and \
                get_item(atom, 'pdbx_PDB_model_num') == '1':
                ndict = {}
                ndict['atom'] = get_item(atom, 'auth_atom_id', conv = int)
                ndict['residue'] = get_item(atom, 'pdbx:auth_comp_id')
                ndict['id'] = int(atom['id'])
                try:
                    ndict['X'] = get_item(atom, 'Cartn_x', conv = float)*scales['X']
                    ndict['Y'] = get_item(atom, 'Cartn_y', conv = float)*scales['Y']
                    ndict['Z'] = get_item(atom, 'Cartn_z', conv = float)*scales['Z']
                except TypeError:
                    continue
                ndict['resnum'] = get_item(atom, 'auth_seq_id', conv = int)
                chain_list.append(deepcopy(ndict))
        
        scheme_list = list()
        for scheme in soup.findAll('pdbx:pdbx_poly_seq_scheme'):
            scheme_list.append({
            'res':get_item(scheme, 'pdbx:auth_mon_id'),
            'authnum':get_item(scheme, 'pdbx:auth_seq_num', conv = int),
            'ndbnum':get_item(scheme, 'pdbx:ndb_seq_num', conv = int),
            'pdbnum':get_item(scheme, 'pdbx:pdb_seq_num', conv = int)
            })
            
        return chain_list, scheme_list

    def calculate_pairwise(self):
        
        pairwise = {}
        chains = list()
        rnum = 0
        for resnum, residues in groupby(self.chain, itemgetter('resnum')):
            reslist = list(residues)
            X = sum([x['X'] for x in reslist])/len(reslist)
            Y = sum([x['Y'] for x in reslist])/len(reslist)
            Z = sum([x['Z'] for x in reslist])/len(reslist)
            chains.append((rnum, X, Y, Z))
            rnum += 1

        for it1, it2 in combinations(chains, 2):
            d = euc_dist(it1[1:], it2[1:])            
            pairwise[(it1[0], it2[0])] = d
            pairwise[(it2[0], it1[0])] = d
            #print it1[0], it2[0]

        self.pairwise = pairwise


    def pairwise_from_seq(self, in_seq):
        if self.pairwise is None:
            self.calculate_pairwise()
        s_seq, i_seq = global_align(self.seq, in_seq)
        scount = -1
        mcount = -1
        i_inds = []
        for s, i in zip(s_seq, i_seq):
            if s != '-':
                scount += 1
            if i != '-':
                mcount += 1
                if scount != -1 and mcount != -1:
                    i_inds.append((mcount, scount))

        pdict = {}
        for s1, s2 in product(i_inds, repeat = 2):
            if s1[1] != s2[1]:
                pdict[(s1[0], s2[0])] = self.pairwise[(s1[1], s2[1])]
        return pdict

def create_scatter(pdb_file, link_file, align_file, out_file, chain):
    
    ref = 'K03455'
    struct_name = gi_from_path(pdb_file)
    prot_name = gi_from_path(align_file)
    
    alignment = Alignment.alignment_from_file(align_file)
    structure = Structure.from_file(pdb_file, chain)
    
    hxb2_seq = alignment.seqs[ref].replace('-', '')
    print 'hxb2', hxb2_seq
    pdict = structure.pairwise_from_seq(hxb2_seq)
    
    ref_nums, align_nums = alignment.convert_numbering(ref)
    
    linkage_dict = {}
    with open(link_file) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            if int(row['Source-End']) - int(row['Source-Start']) == 1 \
                and int(row['Target-End']) - int(row['Target-Start']) == 1 \
                and not row['Correct-Num'].startswith('too'):
                linkage_dict[(int(row['Source-Start']), 
                                int(row['Target-Start']))] = float(row['Total-Score'])
                
    
    fields = ('HXB2Pos1', 'HXB2Pos2', '3dDist', 'Linkage', 
                'Protein', 'Structure')
    with open(out_file, 'w') as handle:
        handle.write('\t'.join(fields)+'\n')
        writer = csv.DictWriter(handle, fieldnames = fields,
                                delimiter = '\t')
        for (p1, p2) in pdict.keys():
            np1 = ref_nums[p1]
            np2 = ref_nums[p2]
            if (np1, np2) in linkage_dict:
                
                writer.writerow({
                'HXB2Pos1':p1,
                'HXB2Pos2':p2,
                '3dDist':pdict[(p1, p2)],
                'Linkage':linkage_dict[(np1, np2)],
                'Protein':prot_name,
                'Structure':struct_name
                })

        
