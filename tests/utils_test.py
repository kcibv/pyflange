import pytest
from pyflange.utils import *
from math import pi
import logging
import pandas as pd



def test_Logger ():
    logger = Logger("pyflange.tests.utils")
    assert isinstance(logger, logging.Logger)
    assert logger.name == "pyflange.tests.utils"
    assert logger.level == LEVEL

    # The same logger is returned for the same id
    assert Logger("pyflange.tests.utils") is logger



def test_load_csv_database ():
    db = load_csv_database("bolts.metric_screws")
    assert isinstance(db, pd.DataFrame)
    assert db["nominal_diameter"]["M16"] == 0.016

    # Results are cached
    assert load_csv_database("bolts.metric_screws") is db

    # Default integer index
    db = load_csv_database("flangesegments.shell_stiffness", index_col=None)
    assert isinstance(db.index, pd.RangeIndex)
    assert "shell_stiffness" in db.columns

    # Unknown database
    with pytest.raises(FileNotFoundError):
        load_csv_database("non_existing_database")



def test_reduce_to_pi ():
    assert reduce_to_pi(0.0) == 0.0
    assert round(reduce_to_pi(0.3), 12) == 0.3
    assert round(reduce_to_pi(-0.3), 12) == 0.3
    assert round(reduce_to_pi(pi), 12) == round(pi, 12)
    assert round(reduce_to_pi(3*pi/2), 12) == round(pi/2, 12)
    assert round(reduce_to_pi(2*pi + 0.3), 12) == 0.3
    assert round(reduce_to_pi(-5*pi/2), 12) == round(pi/2, 12)
