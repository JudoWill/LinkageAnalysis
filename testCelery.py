__author__ = 'will'
from Code.CeleryProcessor import PredictionAnalysis

#redis-server /opt/local/etc/redis.conf


if __name__ == '__main__':
    PredictionAnalysis('HIVData/curated/MergedDir/Env.aln', 'HIVData/curated/MergedDir/Tat.aln', 'tmpres.csv')


