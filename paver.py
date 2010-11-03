from paver.easy import *
import os.path, os
import csv
import ruffus
import urllib2, re
from NCBIUtils import SearchNCBI

options(
    DATA_DIR = 'Data',
    OUT_DIR = 'Results',
)

@task
def touch_data():
    for path, _, files in os.walk(options.DATA_DIR):
        for f in files:
            f = f.replace(' ', '\ ')
            sh('touch %s' % os.path.join(path, f))


@task
def run():

    ruffus.pipeline_run([top_function])

@task
@needs('touch_data', 'run')
def new_run():
    pass


@ruffus.follows('get_sequence_ids')
def top_function():
    pass

@ruffus.posttask(ruffus.touch_file(os.path.join(options.DATA_DIR, 'ListFiles', 'search_sentinal')))
@ruffus.files(os.path.join(options.DATA_DIR, 'ListFiles', 'search_sentinal'),
              os.path.join(options.DATA_DIR, 'ListFiles', 'known_subtypes.list'))
def get_sequence_ids(in_file, out_file):

    search_query = '"Human immunodeficiency virus 1"[porgn] AND 100: 15000[SLEN]'

    id_list = SearchNCBI(search_query)

    with open(out_file, 'w') as handle:
        handle.write('\n'.join(id_list))


