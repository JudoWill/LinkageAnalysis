from paver.easy import *
import os.path, os
import csv
import ruffus

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


@ruffus.files()
def top_function(in_files, out_files):
    pass