"""Internal utilities shared by the other PyFlange modules.

This module collects small helpers that are used throughout the library:

- `Logger`, a logger configured with the PyFlange logging level.
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
