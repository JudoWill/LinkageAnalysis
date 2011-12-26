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

cp LinkageAnalysis/celeryd 	/etc/init.d/
cp LinkageAnalysis/celerydefaults /etc/default/celeryd
chmod +x /etc/init.d/celeryd
/etc/init.d/celeryd start