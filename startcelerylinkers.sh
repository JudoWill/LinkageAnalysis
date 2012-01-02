#!/bin/bash

celeryd --config=clusterceleryconfig --logfile ~/celery.log --loglevel INFO &
disown -a