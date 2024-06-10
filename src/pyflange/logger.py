
import logging
logging.basicConfig()


LEVEL = logging.DEBUG


def Logger (id):
    logger = logging.getLogger(id)
    logger.setLevel(LEVEL)
    return logger


def log_data (object, **kwargs):
    if not hasattr(object, "_data_log"):
        object._data_log = {}
    object._data_log.update(kwargs)


def read_data_log (object, key):
    return object._data_log.get(key) if hasattr(object, "_data_log") else None
