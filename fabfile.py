__author__ = 'will'

from fabric.api import *
from datetime import datetime, timedelta
import time

env.roledefs = {
    'master': ['master'],
    'slaves': ['Node%03i' % x for x in range(1,20)]
}

@roles('master')
def setup_master():
    run('wget http://redis.googlecode.com/files/redis-2.4.5.tar.gz')
    run('tar xzf redis-2.4.5.tar.gz')
    with cd('redis-2.4.5'):
        run('make install')
    with cd('LinkageAnalysis'):
        run('redis-server ./redis.conf')

@roles('master', 'slaves')
def setup_env():
    run('wget http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tgz')
    run('tar xzf Python-2.7.2.tgz')
    with cd('Python-2.7.2'):
        run('./configure')
        run('make install')
    run('wget http://pypi.python.org/packages/source/d/distribute/distribute-0.6.24.tar.gz')
    run('tar xzf distribute-0.6.24.tar.gz')
    with cd('distribute-0.6.24'):
        run('python setup.py install')
    run('wget http://pypi.python.org/packages/source/p/pip/pip-1.0.2.tar.gz')
    run('tar xzf pip-1.0.2.tar.gz')
    with cd('pip-1.0.2'):
        run('python setup.py install')
    with settings(warn_only = True):
        run('git clone git://github.com/JudoWill/LinkageAnalysis.git')
    with cd('LinkageAnalysis'):
        run('pip install -r requirements.pip')

@roles('slaves')
def kill_celery_worker():

    with cd('LinkageAnalysis'):
        run('git reset HEAD --hard')
        run('git pull')
        with settings(warn_only = True):
            lines = run("ps auxww | grep celeryd")
            for line in lines.split('\n'):
                pid = [x for x in line.split() if x.strip()][1]
                run('kill -9 %s' % pid)

@roles('slaves')
def check_workers():
    lastline = run('tail -n 1 celery.log')
    time_format = "%Y-%m-%d %H:%M:%S"
    try:
        tstamp = lastline[1:].split(']')[0].split(',')[0]
        lasttime = datetime.fromtimestamp(time.mktime(time.strptime(tstamp, time_format)))
        tdelta = datetime.now() - lasttime
        print tdelta.total_seconds(), 'seconds since the last process.'
    except IndexError:
        print lastline
    


@roles('master')
def start_celery_workers():

    for slave in env.roledefs['slaves']:
        local("ssh %s 'chmod +x LinkageAnlysis/startcelerylinkers.sh'" % slave)
        local("ssh %s 'LinkageAnlysis/startcelerylinkers.sh'" % slave)
