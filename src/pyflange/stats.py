
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


# References:
# - [1] IEC 61400-6 AMD1 Background document
# - [2] IEC 61400-6:2020 AMD1
'''
This module contains tools for the probabilistic analysis of flanged connections.
In particular, three categories of tools are provided:

- Probability distributions for some flange properties
- Samplers: random realization generators for some flange properties
- Functions for Montecarlo simulations

'''

import scipy.stats as stats
from .bolts import Bolt, Washer, Nut





# =============================================================================
#   
#   PROBABILISTIC DISTRIBUTIONS
#   
# =============================================================================

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





# =============================================================================
#   
#   STANDARD-AGNOSTIC IMPLEMENTATION
#   The function below implement a montecarlo simulation for probabilistic
#   fatigue analysis of flange bolts, in a general form. That is: with no
#   reference to any standard.
#   
# =============================================================================

def is_sampler (value):
    from types import GeneratorType
    return isinstance(value, GeneratorType)


def sampler (random_variable: stats.rv_continuous):
    while True:
        yield random_variable.rvs()


def norm_sampler (mean, cv):
    std = cv * mean
    return sampler( stats.norm(mean, std) )


def lognorm_sampler (mean, coefficient_of_variation):
    from math import log, exp, sqrt
    shape = sqrt( log(coefficient_of_variation**2 + 1) )
    scale = exp(log(mean) - shape**2 / 2)
    return sampler( stats.lognorm(s=shape, loc=0, scale=scale) )


def object_sampler (Class, **params):

    while True:

        params_occurrence = {
            name: next(value) if is_sampler(value) else value
            for name, value in params.items()
        }

        return Class(**params_occurrence)



def bolt_fatigue_sampler (fseg_samp, markov_matrix_samp, fatigue_curve_samp, allowable_damage_samp):
    from .fatigue import BoltFatigueAnalysis

    while True:
        yield BoltFatigueAnalysis(
            fseg = next(fseg_samp),
            flange_mkvm = next(markov_matrix_samp),
            custom_fatigue_curve = next(fatigue_curve_samp),
            allowable_damage = next(allowable_damage_samp) )

# event = {}
#
#         # Extract values from the samples
#         fseg = event["FlangeSegment"] = next(fseg_samp)
#         flange_mkvmat = event["FlangeMarkovMatrix"] = next(markov_matrix_samp)
#         snc = event["FatigueCurve"] = next(fatigue_curve_samp)
#
#         # Transfer the markov matrix to the bolt
#         bending_factor = max(0.5, 0.5 + 0.5*log(fseg.bolt.nominal_diameter/0.036) / log(150/36))
#         bolt_mkvmat = event["BoltMarkovMatrix"] = bolt_markov_matrix(fseg, flange_mkvmat, bending_factor)
#
#         # Generate the next cumulated damage
#         damage = event["Damage"] = snc.cumulated_damage(bolt_mkvmat) * model_uncertainty
#
#         # evaluate the limit state function
#         allowable_damage = event["AllowableDamage"] = next(allowable_damage_samp)
#
#         # evaluate the fatigue life
#         fatigue_life = event["FatigueLife"] = allowable_damage / damage * flange_mkvmat.duration
#
#         yield event


def finite_sampler (sampler, size):
    for i in range(size):
        yield next(sampler)
        

    


# =============================================================================
#   
#   IEC 61400-6 IMPLEMENTATION
#   The function provide a particular montecarlo implementation based on the
#   IEC 61400-6:2025 standard ref. [2] and its background document ref. [1].
#   
# =============================================================================

def standard_gap_size_sampler (flange_diameter, flange_flatness_tolerance):
    from math import pi

    gap_angle_sampler = lognorm_sampler(100/180*pi, 1.0)
    while True:
        gap_angle = next(gap_angle_sampler) % pi
        gap_dist = gap_height_distribution(flange_diameter, flange_flatness_tolerance, gap_angle*flange_diameter/2)
        gap_height = min(gap_dist.rvs(), 2*gap_dist.ppf(0.95))
        yield gap_angle, gap_height



def standard_PolynomialLFlangeSegment_sampler (
        a: float,        # distance between inner face of the flange and center of the bolt hole
        b: float,        # distance between center of the bolt hole and center-line of the shell
        s: float,        # shell thickness
        t: float,        # flange thickness
        R: float,        # shell outer curvature radius
        central_angle: float,     # angle subtended by the flange segment

        Zg: float,       # load applied to the flange segment shell at rest (normally dead weight
                         # of tower + RNA, divided by the number of bolts). Negative if compression.

        bolt: Bolt,      # Bolt object representing the flange segment bolt
        bolt_preload_ratio: float,      # Ratio between mean preload Fp and bolt yield capacity Fy = As*fy.
                                        # According to ref. [1]-section 4:
                                        # - for torque-tightened HV bolts approved for humidity: Fp/Fy = 0.77
                                        # - for torque-tightened HV bolts not approved for humidity: Fp/Fy = 0.62
                                        # - for tension-tightened HV bolts: 0.70 < Fp/Fy < 0.75
                
        bolt_preload_cov: float,        # Coefficient of variation of the bolt preload.
                                        # According to ref. [1]-section 4:
                                        # - for torque-tightened HV bolts approved for humidity: CoV = 0.10
                                        # - for torque-tightened HV bolts not approved for humidity: CoV = 0.15
                                        # - for tension-tightened HV bolts: 0.03 < CoV < 0.05

        Do: float,       # Bolt hole diameter
        washer: Washer,  # Bolt washer
        nut: Nut,        # Bolt nut

        flange_flatness_tolerance: float,

        E: float = 210e9,        # Young modulus of the flange
        G: float = 80.77e9,      # Shear modulus of the flange
        s_ratio: float = 1.0,    # Ratio of bottom shell thickness over s. Default s_botom = s.
        r: float = 0.01,         # Rounding between flange and shell
        k_shell = None           # Custom shell stiffness
    ):
   
    from .flangesegments import PolynomialLFlangeSegment, shell_stiffness
    import numpy as np
    from math import pi
    deg = pi/180

    bolt_yield_capacity = bolt.thread_cross_section.area * bolt.yield_stress
    preload_sampler = norm_sampler(bolt_preload_ratio * bolt_yield_capacity, bolt_preload_cov)
    gap_size_sampler = standard_gap_size_sampler(2*R, flange_flatness_tolerance)
    gap_shape_factor_sampler = norm_sampler(1.0, 0.15)
    tilt_angle_sampler = lognorm_sampler(0.1*deg, 0.50)

    def average_random_preload (preload_sampler, n):
        # Averaging over "50% of the bolts with a gap"
        # I am having a small sub-routine that randomly samples the preload over 
        # half the gap angle and then takes the average of that for each simulation. 
        # Then you get a larger scatter for smaller gap angles.

        n = max(round(n), 1)
        sum = 0
        for i in range(n):
            sum += next(preload_sampler)
        return sum/n

    while True:
        gap_angle, gap_height = next(gap_size_sampler)
        yield PolynomialLFlangeSegment(
                a=a, b=b, s=s, t=t, R=R, central_angle=central_angle, Zg=Zg, bolt=bolt, 
                Do=Do, washer=washer, nut=nut, gap_angle=gap_angle, gap_height=gap_height,
                E=E, G=G, s_ratio=s_ratio, r=r, k_shell=k_shell,

                # Realizations of probabilistic parameters
                Fv = average_random_preload(preload_sampler, gap_angle/central_angle/2),
                gap_shape_factor = next(gap_shape_factor_sampler),   # Factor accounting for a shape different than sinusoidal
                tilt_angle = next(tilt_angle_sampler) % pi           # Flange radia tilt angle
            )



def standard_markov_matrix_sampler (markov_matrix, range_CoV=0.12):
    from .fatigue import MarkovMatrix
    range_coeff_sampler = lognorm_sampler(1, range_CoV)
    while True:
        yield MarkovMatrix(
            cycles = markov_matrix.cycles,
            mean   = markov_matrix.mean,
            range  = markov_matrix.range * next(range_coeff_sampler),
            duration = markov_matrix.duration
        )



def standard_bolt_fatigue_curve_sampler (bolt_nominal_diameter):
    #!!! The fatigue curve sampler must be based on average values, not characteristic
    from .fatigue import BoltFatigueCurve
    stress_factor_samp = norm_sampler(1, 0.10)
    DS_ref_mean = 62e6 # 62 MPa
    while True:
        DS_ref = DS_ref_mean * next(stress_factor_samp)
        yield BoltFatigueCurve(bolt_nominal_diameter, DS_ref, gamma_M=1.0)




