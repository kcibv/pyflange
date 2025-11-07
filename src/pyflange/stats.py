
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
        random_variable (scipy.stats.rv_continuous): any SciPy continuous
            random variable.

    **Returns** a generator that, at every `next` call returns a random realization 
    of the passed random variable.

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
        mean (float): mean value of the normal distribution

        cv (float): coefficient of variation of the normal distribution

    **Returns** a normal distribution sampler with given mean value and
    coefficient of variation.

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
        mean (float): mean value of the log-normal distribution

        cv (float): coefficient of variation of the log-normal distribution

    **Returns** a log-normal distribution sampler with given mean value and
    coefficient of variation.

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


def standard_gap_size_sampler (flange_diameter, flange_flatness_tolerance):
    ''' Sampler that generates random flange gaps according to ref. [1]

    Args:
        flange_diameter (float): the diameter of the flange

        flange_flatness_tolerance (float): the flatness tolerance in mm/mm

    Returns:
        gap_angle (float): random gap angle generated from a log-normal
            sampler with mean value 0.1 deg and CoV = 1.0

        gap_height (float): random gap height, generated from the distribution
            returned by the *gap_height_distribution* function, based on
            the random *gap_angle* and on the given *flange_diameter* and
            *flange_flatness_tolerance* parameters.

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
        k_shell = 'interp'       # Custom shell stiffness
    ):
    '''
    A sampler that yields random pyflange PolynomialLFlangeSegment objects.

    Args:
        a (float): distance between inner face of the flange and center of the 
            bolt hole

        b (float): distance between center of the bolt hole and center-line of 
            the shell

        s (float): shell thickness

        t (float): flange thickness

        R (float): shell outer curvature radius

        central_angle (float): angle subtended by the flange segment

        Zg (float): load applied to the flange segment shell at rest 
            (normally dead weight of tower + RNA, divided by the number of 
            bolts). Negative if compression.

        bolt (pyflange.bolts.Bolt): Bolt object representing the flange 
            segment bolt

        bolt_preload_ratio (float): Ratio between mean preload Fp and bolt 
            yield capacity Fy = As*fy. According to ref. [1]-section 4:
            - for torque-tightened HV bolts approved for humidity: Fp/Fy = 0.77
            - for torque-tightened HV bolts not approved for humidity: Fp/Fy = 0.62
            - for tension-tightened HV bolts: 0.70 < Fp/Fy < 0.75
                
        bolt_preload_cov (float): Coefficient of variation of the bolt preload.
            According to ref. [1]-section 4:
            - for torque-tightened HV bolts approved for humidity: CoV = 0.10
            - for torque-tightened HV bolts not approved for humidity: CoV = 0.15
            - for tension-tightened HV bolts: 0.03 < CoV < 0.05

        Do (float): Bolt hole diameter

        washer (pyflange.bolts.Washer) Bolt washer

        nut (pyflange.bolts.Nut): Bolt nut

        flange_flatness_tolerance (float): flatness tolerance in mm/mm

        E (float): Young modulus of the flange

        G (float): Shear modulus of the flange

        s_ratio (float) Ratio of bottom shell thickness over s. Default s_botom = s.

        r (float): Rounding between flange and shell

        k_shell (str|float|None): Custom shell stiffness. If it equals the string
            'interp' (default), the shell stiffness will be interpolated using
            the funtion *pyflange.flangesemgments.shell_stiffness*; if it is a
            number, the shell stiffness will be taken equal to that number; if it
            is *None*, the shell stiffness will be calculated accordin to the
            simplified formula contained in ref [1].

    *Returns* a sampler that yields random L-Flange segment objects, defined by
    the given deterministic parameters and the following random parameters
    which are generated according to ref. [1].

    - `preload` (`float`), sampled from a Normal distribution with mean value
      given by the passed *bolt_preload_ratio* and CoV given by the passed
      *bolt_preload_cov*.
    - `gap_angle` (`float`) and `gap_height` (`float`), given by the
      *standard_gap_size_generator*, created with flange diameter 2R and
      the given *flange_flatness_tolerance*.
    - `gap_shape_factor` (`float`) sampled from a Normal distribution with
      mean value 1.0 and CoV=0.15
    - `tilt_angle` (`float`) sampled from a log-normal distribution with
      mean value 0.1 deg and CoV=0.50.

    The following example creates a L-Flange segment sampler and generates
    three random flange segments.

    ```py
    from math import pi
    mm = 0.001
    kN = 1000

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

            bolt_preload_ratio = 0.750,     # Ratio between mean preload Fp and Fy = As*fy.
            bolt_preload_cov = 0.03,        # Coefficient of variation of the bolt preload.

            Do = 86*mm,               # Bolt hole diameter
            washer = None,            # Bolt washer
            nut = RoundNut("M80"),    # Bolt nut

            flange_flatness_tolerance = 0.0014,   # 1.4 mm/m

            s_ratio = 1.0       # Ratio of bottom shell thickness over s. Default s_botom = s.
        )

    fseg1 = next(fseg_samp)     # A random L-Flange segment object
    fseg2 = next(fseg_samp)     # Another random L-Flange segment object
    fseg3 = next(fseg_samp)     # Yet another random L-Flange segment object

    ```
    '''
   
    from .flangesegments import PolynomialLFlangeSegment, shell_stiffness
    import numpy as np
    from math import pi
    deg = pi/180

    bolt_yield_capacity = bolt.thread_cross_section.area * bolt.yield_stress
    preload_sampler = norm_sampler(bolt_preload_ratio * bolt_yield_capacity, bolt_preload_cov)
    gap_size_sampler = standard_gap_size_sampler(2*R, flange_flatness_tolerance)
    gap_shape_factor_sampler = norm_sampler(1.0, 0.15)
    tilt_angle_sampler = lognorm_sampler(0.1*deg, 0.50)

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
        gap_angle, gap_height = next(gap_size_sampler)
        yield PolynomialLFlangeSegment(
                a=a, b=b, s=s, t=t, R=R, central_angle=central_angle, Zg=Zg, bolt=bolt, 
                Do=Do, washer=washer, nut=nut, gap_angle=gap_angle, gap_height=gap_height,
                E=E, G=G, s_ratio=s_ratio, r=r, 

                k_shell = np.interp(gap_angle, gap_angles, shell_stiffnesses) if k_shell=='interp' else k_shell,

                # Realizations of probabilistic parameters
                Fv = average_random_preload(preload_sampler, gap_angle/central_angle/2),
                gap_shape_factor = next(gap_shape_factor_sampler),   # Factor accounting for a shape different than sinusoidal
                tilt_angle = next(tilt_angle_sampler) % pi           # Flange radia tilt angle
            )



def standard_markov_matrix_sampler (markov_matrix, range_CoV=0.12):
    ''' Sampler that generates a random markov matrix according to ref [2].

    Args:
        markov_matrix (pyflange.fatigue.MarkovMatrix): the deterministic
            design Markov matrix.

        range_CoV (float): the coefficient of variation of each load
            range, assumed normally-distributed with mean value
            contained in the passed *markov_matrix* parameters.

    **Returns** a MarkovMatrix sampler that yields random Markov matrices
    having the same mean load values and number of cycles as the passed
    *markov_matrix* parameter, but random ranges multiplied by a random
    factor extracted by a normal distribution with mean value 1.0 and
    CoV equal to the passed *range_CoV* parameter.
    '''
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
    ''' Sampler that generates random bolt SN curves, according to ref. [2].

    Args:
        bolt_nominal_diameter (float): the nominal diameter of the bolt

    **Returns** a sampler that yields *pyflange.fatigue.BoltFatigueCurve*
    objects, having a reference stress range randomly generated from a 
    normal distribution with mean value 62 MPa and CoV=0.10.
    '''
    from .fatigue import BoltFatigueCurve
    stress_factor_samp = norm_sampler(1, 0.10)
    DS_ref_mean = 62e6 # 62 MPa
    while True:
        DS_ref = DS_ref_mean * next(stress_factor_samp)
        yield BoltFatigueCurve(bolt_nominal_diameter, DS_ref, gamma_M=1.0)



def fatigue_case_sampler (fseg_samp, markov_matrix_samp, fatigue_curve_samp, allowable_damage_samp):
    ''' Sampler that generates random fatigue cases

    Args:
        fseg_samp (sampler): a *pyflange.flangesegments.FlangeSegment* sampler

        markov_matrix_samp (sampler): a *pyflange.fatigue.MarkovMatric* sampler

        fatigue_curve_samp (sampler): a *pyflange.fatigue.BoltFatigueCurve* sampler

        allowable_damage_samp (sampler): a *float* sampler that generates random allowable damages

    **Returns** a sampler that yields random *pyflange.fatigue.BoltFatigueAnalysis* objects.
    '''
    from .fatigue import BoltFatigueAnalysis

    while True:
        yield BoltFatigueAnalysis(
            fseg = next(fseg_samp),
            flange_mkvm = next(markov_matrix_samp),
            custom_fatigue_curve = next(fatigue_curve_samp),
            allowable_damage = next(allowable_damage_samp) )



