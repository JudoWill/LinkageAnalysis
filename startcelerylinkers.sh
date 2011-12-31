#!/bin/bash

cd LinkageAnalysis
nohup celeryd --config=clusterceleryconfig --logfile ~/celery.log --loglevel INFO &