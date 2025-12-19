
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
'''
This module contains tools for the probabilistic analysis of flanged connections.
In particular, three categories of tools are provided:

- Probability distributions for some flange properties
- Samplers: random realization generators for flange properties

The provided probability distributions are:

- gap_height_distribution

The samplers are just [Python generator functions](https://realpython.com/ref/glossary/generator/)
that yield random values. Samplers can be used in a loop or to generate
random values via the python [next](https://realpython.com/ref/builtin-functions/next/) function.
The following *general-purpose* samplers are available:

- sampler
- norm_sampler
- lognorm_sampler
- fatigue_case_sampler

The following samplers are instead specific to IEC 61400-6:2020:

- standard_gap_size_sampler
- standard_PolynomialLFlangeSegment_sampler
- standard_markov_matrix_sampler
- standard_bolt_fatigue_curve_sampler

An example of how to use the above generator to perform a Montecarlo
simulation, is given in [this example](../examples/montecarlo.md).

The following references are used through this documentation:

- `[1]` IEC 61400-6 AMD1 Background document
- `[2]` IEC 61400-6:2020

'''

from typing import Generator
import scipy.stats as stats
from .bolts import Bolt, Washer, Nut
from math import pi
deg = pi/180





# =============================================================================
#
#   PROBABILISTIC DISTRIBUTIONS
#
# =============================================================================

def gap_height_distribution (flange_diameter, flange_flatness_tolerance, gap_length):
    ''' Evaluates the gap heigh probability distribution according to ref. [1].

    Args:
        flange_diameter (float): The outer diameter of the flange, expressed in meters.
        flange_flatness_tolerance (float): The flatness tolerance, as defined in ref. [1],
            expressed in mm/mm (non-dimensional).
        gap_length (float): The length of the gap, espressed in meters and measured at
            the outer edge of the flange.

    Returns:
        scipy.stats.lognorm: A scipy log-normal variable representing the gap height
            stocastic variable.

    Example:
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
#   SAMPLERS
#   A sampler is just a python generator that yields random values.
#
#   Follows a collection of samplers yielding random realizations of
#   particular flanged connection parameters.
#
# =============================================================================

def sampler (random_variable: stats.rv_continuous):
    ''' Generic distribution-based sampler.

    Args:
        random_variable (scipy.stats.rv_continuous): Any SciPy continuous
            random variable.

    Yields:
        A generator that, at every `next` call returns a random realization
        of the passed random variable.

    Example:
        The following example creates a Normal distribution sampler and generates
        three realizations.

        ```py
        from scipy.stats import norm
        ndist = norm(12.0, 2.0) # normal distribution with mean value 12 and standard deviation 2.

        samp = sampler(ndist)   # sampler based on the ndist distribution

        val1 = next(samp)       # A random value from ndist distribution
        val2 = next(samp)       # Another random value from ndist distribution
        val3 = next(samp)       # Yet another random value from ndist distribution
        ```
    '''
    while True:
        yield random_variable.rvs()


def norm_sampler (mean, cv):
    ''' Sampler based on a Normal distribution.

    Args:
        mean (float): Mean value of the normal distribution.
        cv (float): Coefficient of variation of the normal distribution.

    Returns:
        A normal distribution sampler with given mean value and
        coefficient of variation.

    Example:
        The following example creates a Normal distribution sampler and generates
        three realizations.

        ```py
        samp = norm_sampler(12.0, 0.25)   # normal sampler sampler with mean 12.0
                                          # and CoV = 0.25.

        val1 = next(samp)   # A random value from the normal distribution
        val2 = next(samp)   # Another random value from the normal distribution
        val3 = next(samp)   # Yet another random value from the normal distribution
        ```
    '''
    std = cv * mean
    return sampler( stats.norm(mean, std) )


def lognorm_sampler (mean, cv):
    ''' Sampler based on a Log-Normal distribution.

    Args:
        mean (float): Mean value of the log-normal distribution.
        cv (float): Coefficient of variation of the log-normal distribution.

    Returns:
        A log-normal distribution sampler with given mean value and
        coefficient of variation.

    Example:
        The following example creates a Log-Normal distribution sampler and
        generates three realizations.

        ```py
        samp = lognorm_sampler(12.0, 0.25)  # log-normal sampler sampler with
                                            # mean value 12.0 and CoV = 0.25.

        val1 = next(samp)   # A random value from the log-normal distribution
        val2 = next(samp)   # Another random value from the log-normal distribution
        val3 = next(samp)   # Yet another random value from the log-normal distribution
        ```
    '''
    from math import log, exp, sqrt
    shape = sqrt( log(cv**2 + 1) )
    scale = exp(log(mean) - shape**2 / 2)
    return sampler( stats.lognorm(s=shape, loc=0, scale=scale) )


def standard_gap_sampler (flange_diameter, flange_flatness_tolerance,
                          gap_angle_sampler = lognorm_sampler(100*deg, 1.0),
                          gap_shape_factor_sampler = norm_sampler(1.0, 0.15)):
    ''' Sampler that generates random flange gaps according to ref. [1].

    Args:
        flange_diameter (float): The diameter of the flange.
        flange_flatness_tolerance (float): The flatness tolerance in mm/mm.
        gap_angle_sampler (Generator, optional): A sampler for the gap angle.
            Defaults to a log-normal sampler with mean 0.1 deg and CoV = 1.0.
        gap_shape_factor_sampler (Generator, optional): A sampler for the gap shape factor.
            Defaults to a normal sampler with mean 1.0 and CoV = 0.15.

    Yields:
        tuple: A tuple containing:
            - gap_angle (float): Random gap angle.
            - gap_height (float): Random gap height.

    Example:
        The following example creates a standard gap-size sampler and generate
        three realizations.

        ```py
        # Gap size sampler for 7.5 m falnge with flatness tolerance 1.4 mm/m
        samp = standard_gap_size_sampler(7.5, 0.0014)

        ga1, gh1 = next(samp)  # A random gap-anlge, gap-height pair
        ga2, gh2 = next(samp)  # Another random gap-anlge, gap-height pair
        ga3, gh3 = next(samp)  # Yet another random gap-anlge, gap-height pair

        ```
    '''
    from .flangesegments import Gap

    while True:
        gap_angle = next(gap_angle_sampler) % pi
        gap_dist = gap_height_distribution(flange_diameter, flange_flatness_tolerance, gap_angle*flange_diameter/2)
        gap_height = min(gap_dist.rvs(), 2*gap_dist.ppf(0.95))
        yield Gap(gap_height, gap_angle, next(gap_shape_factor_sampler))



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

        Do: float,       # Bolt hole diameter
        washer: Washer,  # Bolt washer
        nut: Nut,        # Bolt nut

        preload_sampler: Generator,     # Random preload sampler
        gap_sampler: Generator,         # gap object random sampler

        tilt_sampler: Generator = lognorm_sampler(0.1*deg, 0.50),

        E: float = 210e9,        # Young modulus of the flange
        G: float = 80.77e9,      # Shear modulus of the flange
        s_ratio: float = 1.0,    # Ratio of bottom shell thickness over s. Default s_botom = s.
        r: float = 0.01,         # Rounding between flange and shell
        k_shell = 'interp'       # Custom shell stiffness
    ):
    ''' A sampler that yields random pyflange PolynomialLFlangeSegment objects.

    Args:
        a (float): Distance between inner face of the flange and center of the
            bolt hole.
        b (float): Distance between center of the bolt hole and center-line of
            the shell.
        s (float): Shell thickness.
        t (float): Flange thickness.
        R (float): Shell outer curvature radius.
        central_angle (float): Angle subtended by the flange segment.
        Zg (float): Load applied to the flange segment shell at rest
            (normally dead weight of tower + RNA, divided by the number of
            bolts). Negative if compression.
        bolt (pyflange.bolts.Bolt): Bolt object representing the flange
            segment bolt.
        Do (float): Bolt hole diameter.
        washer (pyflange.bolts.Washer): Bolt washer.
        nut (pyflange.bolts.Nut): Bolt nut.
        preload_sampler (Generator): Random preload sampler.
        gap_sampler (Generator): Gap object random sampler.
        tilt_sampler (Generator, optional): Random tilt angle sampler.
            Defaults to lognorm_sampler(0.1*deg, 0.50).
        E (float, optional): Young modulus of the flange. Defaults to 210e9.
        G (float, optional): Shear modulus of the flange. Defaults to 80.77e9.
        s_ratio (float, optional): Ratio of bottom shell thickness over s.
            Defaults to 1.0.
        r (float, optional): Rounding between flange and shell. Defaults to 0.01.
        k_shell (str|float|None, optional): Custom shell stiffness.
            If 'interp' (default), stiffness is interpolated.
            If a number, it's used as the stiffness.
            If None, a simplified formula is used.

    Yields:
        pyflange.flangesegments.PolynomialLFlangeSegment: A random L-Flange segment object.

    Example:
        The following example creates a L-Flange segment sampler and generates
        three random flange segments.

        ```py
        from math import pi
        mm = 0.001
        kN = 1000

        import pyflange.stats as stats
        from pyflange.bolts import StandardMetricBolt, RoundNut
        N_BOLTS = 156

        fseg_samp = standard_PolynomialLFlangeSegment_sampler (
                a = 150*mm,        # distance between inner face of the flange and center of the bolt hole
                b = 122*mm,        # distance between center of the bolt hole and center-line of the shell
                s =  54*mm,        # shell thickness
                t = 172*mm,        # flange thickness
                R = 4000*mm,       # shell outer curvature radius
                central_angle = 2*pi / N_BOLTS,  # angle subtended by the flange segment

                Zg = -18044*kN / N_BOLTS,      # load applied to the flange segment shell at rest

                # Bolt object representing the flange segment bolt
                bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm,
                        shank_diameter_ratio=76.1/80, stud=True),

                Do = 86*mm,               # Bolt hole diameter
                washer = None,            # Bolt washer
                nut = RoundNut("M80"),    # Bolt nut

                # bolt preload random sampler
                preload_sampler = stats.norm_sampler(2932.24*kN, 0.03), 

                # gap object random sampler
                gap_sampler = stats.standard_gap_sampler(
                            flange_diameter = 8000*mm,
                            flange_flatness_tolerance = 0.0014,   # 1.4 mm/m
                            gap_angle_sampler = stats.lognorm_sampler(100*deg, 1.0),
                            gap_shape_factor_sampler = stats.norm_sampler(1.0, 0.15)
                        ),

                s_ratio = 1.0       # Ratio of bottom shell thickness over s. Default s_botom = s.
            )

        fseg1 = next(fseg_samp)     # A random L-Flange segment object
        fseg2 = next(fseg_samp)     # Another random L-Flange segment object
        fseg3 = next(fseg_samp)     # Yet another random L-Flange segment object

        ```
    '''

    from .flangesegments import PolynomialLFlangeSegment, shell_stiffness
    import numpy as np

    # stiffness interpolaion
    gap_angles = np.linspace(10*deg, 180*deg, 100)
    shell_stiffnesses = np.array([shell_stiffness(R, s, gap_angle) for gap_angle in gap_angles])

    # preload averaging over gap
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

    # generate random flange segments
    while True:
        gap = next(gap_sampler)
        yield PolynomialLFlangeSegment(
                a=a, b=b, s=s, t=t, R=R, central_angle=central_angle, Zg=Zg, bolt=bolt,
                Do=Do, washer=washer, nut=nut, gap=gap, E=E, G=G, s_ratio=s_ratio, r=r,

                k_shell = np.interp(gap.angle, gap_angles, shell_stiffnesses) if k_shell=='interp' else k_shell,

                # Realizations of probabilistic parameters
                Fv = average_random_preload(preload_sampler, gap.angle/central_angle/2),
                tilt_angle = next(tilt_sampler) % pi           # Flange radia tilt angle
            )



def standard_markov_matrix_sampler (markov_matrix, mean_range_coeff=1.0, range_CoV=0.12):
    ''' Sampler that generates a random markov matrix according to ref [2].

    Args:
        markov_matrix (pyflange.fatigue.MarkovMatrix): The deterministic
            design Markov matrix.
        mean_range_coeff (float, optional): The mean coefficient of each load
            range, assumed log-normally-distributed with mean value
            contained in the passed *markov_matrix* parameters. Defaults to 1.0.
        range_CoV (float, optional): The coefficient of variation of each load
            range, assumed log-normally-distributed with mean value
            contained in the passed *markov_matrix* parameters. Defaults to 0.12.

    Yields:
        pyflange.fatigue.MarkovMatrix: A random Markov matrix.
    '''
    from .fatigue import MarkovMatrix
    range_coeff_sampler = lognorm_sampler(mean_range_coeff, range_CoV)
    while True:
        yield MarkovMatrix(
            cycles = markov_matrix.cycles,
            mean   = markov_matrix.mean,
            range  = markov_matrix.range * next(range_coeff_sampler),
            duration = markov_matrix.duration
        )



def standard_bolt_fatigue_curve_sampler (bolt_nominal_diameter, mean_stress_factor=1.0, stress_factor_CoV=0.10):
    ''' Sampler that generates random bolt SN curves, according to ref. [2].

    Args:
        bolt_nominal_diameter (float): The nominal diameter of the bolt.
        mean_stress_factor (float, optional): The mean value of the fatigue class
            coefficinet, assumed normally-distributed. Defaults to 1.
        range_CoV (float, optional): The coefficient of variation of the fatigue
            class, assumed normally distributed. It defaults to 0.10.

    Yields:
        pyflange.fatigue.BoltFatigueCurve: A random bolt fatigue curve object.
    '''
    from .fatigue import BoltFatigueCurve
    stress_factor_samp = norm_sampler(mean_stress_factor, stress_factor_CoV)
    DS_ref_mean = 62e6 # 62 MPa
    while True:
        DS_ref = DS_ref_mean * next(stress_factor_samp)
        yield BoltFatigueCurve(bolt_nominal_diameter, DS_ref, gamma_M=1.0)



def fatigue_case_sampler (fseg_samp, markov_matrix_samp, fatigue_curve_samp, allowable_damage_samp):
    ''' Sampler that generates random fatigue cases.

    Args:
        fseg_samp (Generator): A `pyflange.flangesegments.FlangeSegment` sampler.
        markov_matrix_samp (Generator): A `pyflange.fatigue.MarkovMatrix` sampler.
        fatigue_curve_samp (Generator): A `pyflange.fatigue.BoltFatigueCurve` sampler.
        allowable_damage_samp (Generator): A `float` sampler that generates random allowable damages.

    Yields:
        pyflange.fatigue.BoltFatigueAnalysis: A random `BoltFatigueAnalysis` object.
    '''
    from .fatigue import BoltFatigueAnalysis

    while True:
        yield BoltFatigueAnalysis(
            fseg = next(fseg_samp),
            flange_mkvm = next(markov_matrix_samp),
            custom_fatigue_curve = next(fatigue_curve_samp),
            allowable_damage = next(allowable_damage_samp) )
