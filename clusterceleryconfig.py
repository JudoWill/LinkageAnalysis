BROKER_BACKEND = 'redis'
BROKER_URL = 'redis://master:6379'
CELERY_RESULT_BACKEND = "redis"
CELERY_REDIS_HOST = "master"
CELERY_REDIS_PORT = 6379
CELERY_REDIS_DB = 0
CELERY_DISABLE_RATE_LIMITS = True
CELERY_IMPORTS = ('Code.CeleryProcessor',)
CELERYD_TASK_TIME_LIMIT = 60*7
CELERYD_SOFT_TASK_TIME_LIMIT = 20
CELERYD_PREFETCH_MULTIPLIER = 1
CELERYD_CONCURRENCY = 20