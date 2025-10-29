import pytest
import pyflange.stats as stats
from math import pi


def test_gap_height_distribution ():
    D = 7.5
    u = 0.0014
    L = pi/6 * D/2

    gap_dist = stats.gap_height_distribution(D, u, L)

    assert round(gap_dist.mean(), 6) == 0.000288
    assert round(gap_dist.std(), 6) == 0.000350
    assert round(gap_dist.ppf(0.95), 6) == 0.000876
