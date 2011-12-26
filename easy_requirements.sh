#!/bin/bash
easy_install paver
easy_install ruffus
easy_install nose
easy_install memorised
easy_install python-memcached
easy_install pyyaml
easy_install BeautifulSoup
easy_install suds
easy_install dendropy
easy_install pylru
easy_install redis
easy_install celery

cd LinkageAnalysys
nohup celeryd --autoscale=10,3 --config=clusterceleryconfig &
