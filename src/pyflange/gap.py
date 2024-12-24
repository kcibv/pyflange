# pyFlange - python library for large flanges design
# Copyright (C) 2024  KCI The Engineers B.V.,
#                     Siemens Gamesa Renewable Energy B.V.,
#                     Nederlandse Organisatie voor toegepast-natuurwetenschappelijk onderzoek TNO.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License version 3 for more details.
#
# You should have received a copy of the GNU General Public License
# version 3 along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
The ``gap`` module contains tools for modelling flange gaps.

It corrently contains only one gap model, that is the sinusoidal gap, as defined
in ref. [1], section 6.7.5.2.


**REFERENCES**

- **[1]** IEC 61400-6:2020/AMD1:2024 - Wind Energy Generation Systems - Part 6: Tower and foundation design requirements
'''


def gap_height_distribution (flange_diameter, flange_flatness_tolerance, gap_length):
    ''' Evaluate the gap heigh probability distribution according to ref. [1].

    Args:
        flange_diameter (float): The outer diameter of the flange, expressed in meters.

        flange_flatness_tolerance (float): The flatness tolerance, as defined in ref. [1],
            expressed in mm/mm (non-dimensional).

        gap_length (float): The length of the gap, espressed in meters and measured at
            the outer edge of the flange.

    Returns:
        dist (scipy.stats.lognorm): a [scipy log-normal variable](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.lognorm.html),
            representing the gap height stocastic variable.

    The following example, creates a gap distribution and the calculates the 95% quantile
    of the gap height

    ```python
    from pyflange.gap import gap_height_distribution

    D = 7.50      # Flange diameter in meters
    u = 0.0014    # Flatness tolerance (non-dimensional)
    L = 1.22      # Gap length
    gap_dist = gap_height_distribution(D, u, L)     # a lognorm distribution object

    u95 = gap_dist.ppf(0.95)    # PPF is the inverse of CDF. See scipy.stats.lognorm documentation.
    ```
    '''

    from math import pi, log, exp, sqrt
    from scipy.stats import lognorm

    k_mean = (6.5/flange_diameter * (flange_flatness_tolerance/0.0014) * (0.025*gap_length**2 + 0.12*gap_length)) / 1000
    gap_angle_deg = (gap_length / (flange_diameter/2)) / pi*180
    k_COV = 0.35 + 200 * gap_angle_deg**(-1.6)
    k_std = k_mean * k_COV

    shape = sqrt( log(k_COV**2 + 1) )
    scale = exp(log(k_mean) - shape**2 / 2)

    return lognorm(s=shape, loc=0, scale=scale)
