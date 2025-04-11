

# ------------------------------------------------------------------------------
# Logging utilities
# ------------------------------------------------------------------------------

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




# ------------------------------------------------------------------------------
#  CSV Databases
# ------------------------------------------------------------------------------
from functools import cache

@cache
def load_csv_database (db_name, index_col="id"):
    import os
    import pandas as pd
    db_path = os.path.join(os.path.dirname(__file__), f"data/{db_name}.csv")
    return pd.read_csv(db_path, index_col=index_col)
