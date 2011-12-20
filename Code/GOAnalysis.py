from __future__ import division
import csv, os.path, os
from collections import defaultdict
from itertools import combinations, imap, product, ifilter, repeat
from operator import methodcaller, contains
from functools import partial
from copy import deepcopy

def GO2name(term):
    try:    
        if '~' in term:
            return term.split('~')[1]
        else:
            return term.split(':')[1]
    except IndexError:
        return term.strip()


def transposedict(in_d):
    res = defaultdict(set)
    for key, vals in in_d.items():
        for val in vals:
            res[val].add(key)
    return res

def calculate_term_overlap(term_dict, old_dict = None, check_terms = None):
    if old_dict is not None:
        vals = old_dict
    else:
        vals = defaultdict(int)

    if check_terms is None:
        it = combinations(term_dict.keys(), 2)
    else:
        it = product(term_dict.keys(), check_terms)

    for k1, k2 in it:
        inter = term_dict[k1] & term_dict[k2]
        un = term_dict[k1] | term_dict[k2]
        vals[(k1, k2)] = (len(inter)/len(un), len(un))
    
    return vals

def check_subsets(term_dict):
    vals = {}
    for k1, k2 in product(term_dict.keys()):
        vals[(k1, k2)] = vals[k2] < vals[k1]
    return vals

def score_mapping(mapping, gene2term):
    
    overlap_count = 0
    for i1, i2 in combinations(mapping.values(), 2):
        overlap_count += len(i1 & i2)

    all_terms = reduce(lambda x,y: x|y, mapping.values(), set())
    unmapped_count = len(set(gene2term.keys()) - all_terms)

    return (overlap_count + unmapped_count)/len(gene2term)
        

def join_terms(terms, term2gene, term_overlap):

    def split_terms(terms):
        for term in terms:
            for it in term.split('||'):
                yield it
                    
    nterm = '||'.join(sorted(list(split_terms(terms))))
    pop_terms = []

    for term in terms:
        term2gene[nterm] |= term2gene.pop(term)
        pop_terms += [x for x in term_overlap.iterkeys() if term in x]

    for term in pop_terms:
        term_overlap.pop(term, None)
 
    return term2gene, calculate_term_overlap(term2gene, old_dict = term_overlap, 
                                                check_terms = (nterm,))



def make_mapping(term2gene, term_overlap, num_left, num_check, total_items = 15):
    
    def score_fun(arg):
        key, item = arg        
        key_join = key[0]+key[1]
        num_terms = sum([k == '|' for k in key_join])/2
        return num_terms, -item[1], -item[0], key


    if num_left == 1 or num_check == 1:
        mapping = term2gene
        while len(mapping) > total_items:
            _, _, _, (k1, k2) = min(imap(score_fun, term_overlap.iteritems()))
            print len(mapping), 'joining', k1, 'WITH', k2
            mapping, term_overlap = join_terms((k1,k2), mapping, term_overlap)
        score = score_mapping(mapping, transposedict(mapping))
        return mapping, score
       

    items = sorted(imap(score_fun, term_overlap.iteritems()))
    mscore = None
    for _, _, _, (k1, k2) in items[:num_check]:
        print 'joining', k1, 'WITH', k2
        jmapping, joverlap = join_terms((k1,k2), term2gene, term_overlap)
        nmapping, nscore = make_mapping(jmapping, joverlap, num_left-1, max(num_check-1,1))
        print num_left, nscore, k1, k2
        if mscore is None or nscore < mscore:
            mmapping = deepcopy(nmapping)
            mscore = nscore

    return mmapping
            
    
    




if __name__ == '__main__':
    

    num_terms = 15
    terms = ('GOTERM_MF_FAT', 'GOTERM_BP_FAT', 'KEGG_PATHWAY')
    for term in terms:
        gene2term = defaultdict(set)

        with open('GOTerms.txt') as handle:
            for row in csv.DictReader(handle, delimiter = '\t'):
                if row[term]:
                    gene2term[row['ID']] |= set([GO2name(x) for x in row[term].split(',') if x])

        term2gene = transposedict(gene2term)
        term_overlap = calculate_term_overlap(term2gene)
        
        mapping, score = make_mapping(term2gene, term_overlap, 15, 1)
        print score, mapping.keys()
    
