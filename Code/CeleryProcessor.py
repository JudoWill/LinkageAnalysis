from __future__ import division
from celery.task import task
import LinkUtils
import LinkFields
import AlignUtils
from AlignUtils import Alignment, prot_from_path
from itertools import product, izip, count, combinations, islice
from Queue import Queue
from threading import Thread
from functools import partial
import csv


@task()
def link_calculator(row, submats, seq1, seq2):

    c1 = AlignUtils.make_counts(seq1)
    c2 = AlignUtils.make_counts(seq1)

    row['S1-Entropy'] = AlignUtils.calculate_entropy(seq1)
    row['S2-Entropy'] = AlignUtils.calculate_entropy(seq2)
    row['S12-Mapping'] = LinkUtils.prediction_mapping(seq1, seq2)
    row['S21-Mapping'] = LinkUtils.prediction_mapping(seq2, seq1)
    row['SeqLength'] = len(seq1)
    row['S1-Cons'] = max(x/len(seq1) for x in c1.values())
    row['S2-Cons'] = max(x/len(seq2) for x in c2.values())

    processfuns = []
    for name, mat in submats:
        processfuns.append(('SBASC_'+name,
            partial(LinkUtils.calculate_SBASC, mat)))

    processfuns.append(('Mutual_Info', LinkUtils.calculate_mutual_info))
    processfuns.append(('OMES', LinkUtils.calculate_OMES))
    suffs = ['_raw', '_pval', '_null', '_count']

    for name, func in processfuns:
        res = LinkUtils.calculate_vals(seq1, seq2, func, minreps=20, maxreps = 100)
        for val, suff in zip(res, suffs):
            row[name+suff] = val

    return row


def task_loader(que, a1, a2, defaults, submats, minseqs, issame):

    calc = LinkUtils.LinkCalculator()
    headers = sorted(set(a1.seqs.keys()) & set(a2.seqs.keys()))

    s1cols = list(a1.iterate_columns(headers))
    s2cols = list(a2.iterate_columns(headers))

    if issame:
        iterable = combinations(range(len(s2cols)), 2)
    else:
        iterable = product(range(len(s1cols)), range(len(s2cols)))
    for ind1, ind2 in iterable:
        #print 'trying', ind1, ind2
        row = {}
        row.update(defaults)
        try:
            seq1 = s1cols[ind1]
            seq2 = s2cols[ind2]
            row['S1-Start'] = ind1
            row['S1-End'] = ind1+1
            row['S2-Start'] = ind2
            row['S2-End'] = ind2+1
        except IndexError:
            continue
        cseq1 = ''
        cseq2 = ''
        for s1, s2 in zip(seq1, seq2):
            if s1 != '-' and s2 != '-':
                cseq1 += s1
                cseq2 += s2

        if len(cseq1) > minseqs:
            #print 'loaded', ind1, ind2
            que.put(link_calculator.delay(row, submats, cseq1, cseq2))
    que.put(None)


def convert_row_to_writeable_rows(row, rmfields):

    s12row = {}
    s21row = {}
    s12row.update(row)
    s21row.update(row)
    for key in row.keys():
        if key.startswith('S1-'):
            s12row['Source-'+key[3:]] = row[key]
            s21row['Target-'+key[3:]] = row[key]
        elif key.startswith('S2-'):
            s12row['Target-'+key[3:]] = row[key]
            s21row['Source-'+key[3:]] = row[key]

    s12row['Total-Score'] = sum(z for _, _, z in row['S12-Mapping'])/row['SeqLength']
    s21row['Total-Score'] = sum(z for _, _, z in row['S21-Mapping'])/row['SeqLength']
    for source, target, val in s12row['S12-Mapping']:
        s12row.update({'Source-Seq':source,
                    'Target-Seq':target,
                    'Correct-Num':val})
        yield s12row
        s12row.update(rmfields)

    for source, target, val in s21row['S21-Mapping']:
        s21row.update({'Source-Seq':source,
                       'Target-Seq':target,
                       'Correct-Num':val})
        yield s21row
        s21row.update(rmfields)




def PredictionAnalysis(align1, align2, outfile, cons_cut = 0.99, **kwargs):

    a1 = Alignment.alignment_from_file(align1)
    a2 = Alignment.alignment_from_file(align2)

    print 'loaded alignments'
    sprot = prot_from_path(align1)
    tprot = prot_from_path(align2)

    defaults = dict(zip(LinkFields.LINK_FIELDS, [None]*len(LinkFields.LINK_FIELDS)))
    defaults['S1-Prot']=sprot
    defaults['S2-Prot']=tprot
    defaults.pop('Source-Start')
    defaults.pop('Source-End')
    defaults.pop('Target-Start')
    defaults.pop('Target-End')

    calculator = LinkUtils.LinkCalculator()
    rmheaders = dict((head, None) for head in calculator.get_fields())

    submats = LinkUtils.get_all_sub_mats()

    process_que = Queue(1000)
    loader = Thread(target=task_loader,
                    args= (process_que, a1, a2, defaults,submats, 50, align1==align2))
    loader.start()

    ohandle = open(outfile, 'w')
    owriter = csv.DictWriter(ohandle, LinkFields.LINK_FIELDS,
        delimiter = '\t', extrasaction='ignore')
    owriter.writerow(dict(zip(LinkFields.LINK_FIELDS,
        LinkFields.LINK_FIELDS)))
    print 'waiting for first'
    item = process_que.get()
    while item is not None:
        row = item.get()
        print row['S1-Start'], row['S2-End']
        owriter.writerows(convert_row_to_writeable_rows(row, rmheaders))
        item = process_que.get()





