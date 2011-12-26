#!/bin/bash

wget http://www.python.org/ftp/python/2.7.2/Python-2.7.2.tgz
tar xzf Python-2.7.2.tgz
cd Python-2.7.2
./configure
make install
cd ~/

wget http://pypi.python.org/packages/source/d/distribute/distribute-0.6.24.tar.gz
tar xzf distribute-0.6.24.tar.gz
cd distribute-0.6.24
python setup.py install
cd ~/

wget http://pypi.python.org/packages/source/p/pip/pip-1.0.2.tar.gz
tar xzf pip-1.0.2.tar.gz
cd pip-1.0.2
python setup.py install
cd ~/

wget http://redis.googlecode.com/files/redis-2.4.5.tar.gz
tar xzf redis-2.4.5.tar.gz
cd redis-2.4.5
make
cd ~/

cd LinkageAnalysis
redis-server redis.conf
pip install -r requirements.pip

nohup celeryd --config=cluterceleryconfig --autoscale=10,3 --logfile ~/celery.log --loglevel INFO &

