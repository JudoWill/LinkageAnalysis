#!/bin/bash

cd LinkageAnalysis
celeryd --config=clusterceleryconfig --logfile ~/celery.log --loglevel INFO &
disown -a