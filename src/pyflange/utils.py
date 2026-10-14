"""Internal utilities shared by the other PyFlange modules.

This module collects small helpers that are used throughout the library:

- Logging helpers: `Logger`, `log_data` and `read_data_log`.
- `load_csv_database`, which loads the tabular data (bolt dimensions, material
  properties, shell stiffness results, etc.) shipped with the package.
- `reduce_to_pi`, an angle reduction function.

These helpers are mainly meant for internal use, but they are public and
documented because they may be useful when extending the library.
"""

# ------------------------------------------------------------------------------
# Logging utilities
# ------------------------------------------------------------------------------

import logging
logging.basicConfig()


LEVEL = logging.DEBUG


def Logger (id):
    """Create a logger configured with the PyFlange logging level.

    Args:
        id (str): Name of the logger. By convention, the `__name__` of the
            calling module.

    Returns:
        logging.Logger: The logger registered under the given name, with its
            level set to the module-level `LEVEL` constant (`logging.DEBUG`
            by default).

    Example:
        ```python
        from pyflange.utils import Logger

        logger = Logger(__name__)
        logger.warning("Something looks suspicious")
        ```
    """
    logger = logging.getLogger(id)
    logger.setLevel(LEVEL)
    return logger


def log_data (object, **kwargs):
    """Attach intermediate calculation values to an object.

    The values are stored in a private dictionary on the object (created on
    first use), so that intermediate results of a calculation can be inspected
    afterwards with `read_data_log`. Existing keys are overwritten.

    Args:
        object (Any): The object the data is attached to.
        **kwargs: Named values to log.

    Example:
        ```python
        from pyflange.utils import log_data, read_data_log

        class Calc: pass
        calc = Calc()
        log_data(calc, k_shell=1.5e9, f_tot=2.0)
        read_data_log(calc, "f_tot")    # 2.0
        ```
    """
    if not hasattr(object, "_data_log"):
        object._data_log = {}
    object._data_log.update(kwargs)


def read_data_log (object, key):
    """Read a value previously stored with `log_data`.

    Args:
        object (Any): The object the data was attached to.
        key (str): Name of the logged value.

    Returns:
        Any: The logged value, or `None` if nothing has been logged on the
            object or if the key is unknown.
    """
    return object._data_log.get(key) if hasattr(object, "_data_log") else None




# ------------------------------------------------------------------------------
#  CSV Databases
# ------------------------------------------------------------------------------
from functools import cache

@cache
def load_csv_database (db_name, index_col="id"):
    """Load one of the CSV databases shipped with the package.

    The files are read from the `data` folder of the package, as
    `data/<db_name>.csv`. Results are cached, so each database is read from
    disk only once.

    Args:
        db_name (str): Name of the database, without the `.csv` extension
            (e.g. `"bolts.metric_screws"`).
        index_col (str | None, optional): Name of the column used as row index.
            Use `None` to keep the default integer index. Defaults to `"id"`.

    Returns:
        pandas.DataFrame: The database content. Since the result is cached,
            it should be treated as read-only.

    Raises:
        FileNotFoundError: If no database with the given name exists.
    """
    import os
    import pandas as pd
    db_path = os.path.join(os.path.dirname(__file__), f"data/{db_name}.csv")
    return pd.read_csv(db_path, index_col=index_col)



# ------------------------------------------------------------------------------
#  Angle Reduction
# ------------------------------------------------------------------------------

def reduce_to_pi (angle):
    """Reduce an angle to the interval [0, pi].

    The angle is first wrapped to (-pi, pi] and then its absolute value is
    taken. This is useful for angles whose sign is irrelevant, such as a tilt
    angle magnitude.

    Args:
        angle (float): The angle, in radians.

    Returns:
        float: The equivalent non-negative angle in radians, between 0 and pi.

    Example:
        ```python
        from math import pi
        reduce_to_pi(3*pi/2)    # pi/2
        reduce_to_pi(-0.3)      # 0.3
        ```
    """
    from math import sin, cos, atan2
    return abs(atan2(sin(angle), cos(angle)))
