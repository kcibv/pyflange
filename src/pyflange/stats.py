
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

'''Probability distributions and random samplers for probabilistic flange analyses.

This module contains tools for the probabilistic analysis of flanged
connections. In particular, two categories of tools are provided:

- Probability distributions for some flange properties.
- Samplers: random realization generators for flange properties.

The provided probability distributions are:

- `gap_height_distribution` (IEC 61400-6 AMD1)
- `gap_height_distribution_ED2` and `gap_angle_distribution_ED2` (IEC 61400-6 ED2)

The samplers are just [Python generator functions](https://realpython.com/ref/glossary/generator/)
that yield random values. Samplers can be used in a loop or to generate
random values via the python [next](https://realpython.com/ref/builtin-functions/next/) function.
The following *general-purpose* samplers are available:

- `sampler`
- `norm_sampler`
- `lognorm_sampler`
- `weibull_sampler`
- `beta_sampler`
- `product_sampler`
- `clayton_copula_sampler`
- `clayton_gap_sampler`
- `fatigue_case_sampler`

The following samplers are instead specific to IEC 61400-6:2020:

- `standard_gap_sampler`
- `standard_gap_sampler_ED2`
- `standard_PolynomialLFlangeSegment_sampler`
- `standard_markov_matrix_sampler`
- `standard_bolt_fatigue_curve_sampler`

An example of how to use the above generators to perform a Montecarlo
simulation, is given in [this example](../examples/montecarlo.md).

The following references are used through this documentation:

- `[1]` IEC 61400-6 AMD1 Background document
- `[2]` IEC 61400-6:2020
- `[3]` IEC 61400-6 ED2 (draft)

'''

from typing import Generator
import scipy.stats as stats
from .bolts import Bolt, Washer, Nut
from metrum.units import deg, MPa





# =============================================================================
#
#   PROBABILISTIC DISTRIBUTIONS
#
# =============================================================================

def gap_height_distribution (flange_diameter, flange_flatness_tolerance, gap_length):
    ''' Evaluates the gap height probability distribution according to ref. [1].

    The gap height is modeled as a log-normal random variable, whose mean value
    and coefficient of variation depend on the flange diameter, the flatness
    tolerance and the gap length.

    Args:
        flange_diameter (float): The outer diameter of the flange, expressed in meters.
        flange_flatness_tolerance (float): The flatness tolerance, as defined in ref. [1],
            expressed in mm/mm (non-dimensional).
        gap_length (float): The length of the gap, expressed in meters and measured at
            the outer edge of the flange.

    Returns:
        scipy.stats.rv_frozen: A frozen SciPy log-normal distribution
            (`scipy.stats.lognorm`) representing the gap height stochastic
            variable, expressed in meters.

    Example:
        The following example creates a gap distribution and then calculates
        the 95% quantile of the gap height.

        ```python
        from pyflange.stats import gap_height_distribution

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



def gap_height_distribution_ED2 (flange_diameter, flange_flatness_tolerance, gap_length):
    ''' Evaluates the gap height probability distribution according to ref. [3].

    The gap height is modeled as a Weibull random variable, whose shape and
    scale parameters depend on the flange diameter, the flatness tolerance and
    the gap angle.

    Args:
        flange_diameter (float): The outer diameter of the flange, in meters.
        flange_flatness_tolerance (float): The flatness tolerance, in mm/mm
            (non-dimensional).
        gap_length (float): The length of the gap, in meters, measured at the
            outer edge of the flange.

    Returns:
        scipy.stats.rv_frozen: A frozen SciPy Weibull distribution
            (`scipy.stats.weibull_min`) representing the gap height, in meters.

    Example:
        ```python
        from math import pi
        from pyflange.stats import gap_height_distribution_ED2

        D = 7.50                         # Flange diameter in meters
        u = 0.0014                       # Flatness tolerance (1.4 mm/m)
        L = pi/6 * D/2                   # Length of a 30 deg gap
        dist = gap_height_distribution_ED2(D, u, L)
        u95 = dist.ppf(0.95)             # 95% quantile of the gap height
        ```
    '''
    from math import pi
    from scipy.stats import weibull_min

    # Convert parameters to units used in ED2 formulae
    u_tol = flange_flatness_tolerance / 0.001             # mm/m
    gap_angle = gap_length / (flange_diameter/2) / deg    # degrees

    # Weibull shape factor
    alpha = 0.667 + 0.0232 * gap_angle - 0.0000489 * gap_angle**2

    # Weibull scale factor, in mm
    beta1 = (0.00588 + 0.000304*flange_diameter) * u_tol
    beta2 = (0.0000109 * u_tol**2 + 0.0000378 * u_tol) * (0.2274*flange_diameter) \
            * (1.0447 * (pi/12 * flange_diameter * u_tol)**(-0.291))
    beta = beta1 * gap_angle + beta2 * gap_angle**2

    return weibull_min(alpha, scale=beta*0.001)


def gap_angle_distribution_ED2 (flange_diameter, flange_flatness_tolerance):
    ''' Evaluates the gap angle probability distribution according to ref. [3].

    The gap angle is modeled as a Beta random variable defined over the
    interval [0, pi], whose parameters depend on the flange diameter and
    flatness tolerance.

    Args:
        flange_diameter (float): The outer diameter of the flange, in meters.
        flange_flatness_tolerance (float): The flatness tolerance, in mm/mm
            (non-dimensional).

    Returns:
        scipy.stats.rv_frozen: A frozen SciPy Beta distribution
            (`scipy.stats.beta`) representing the gap angle, in radians.
    '''
    from math import pi
    from scipy.stats import beta

    u_tol = flange_flatness_tolerance / 0.001     # mm/m

    a = 0.9 * 1.103 * (pi * flange_diameter * u_tol/12)**(-0.216)
    b = 1.3373 * (flange_diameter * u_tol**2)**(-0.146)

    return beta(a, b, loc=0, scale=pi)





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
        float: A random realization of the passed random variable, at every
            `next` call.

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
        Generator: A normal distribution sampler with given mean value and
            coefficient of variation.

    Example:
        The following example creates a Normal distribution sampler and generates
        three realizations.

        ```py
        samp = norm_sampler(12.0, 0.25)   # normal sampler with mean 12.0
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
        Generator: A log-normal distribution sampler with given mean value
            and coefficient of variation.

    Example:
        The following example creates a Log-Normal distribution sampler and
        generates three realizations.

        ```py
        samp = lognorm_sampler(12.0, 0.25)  # log-normal sampler with
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


def weibull_sampler (shape, scale=1.0):
    ''' Sampler based on a Weibull distribution.

    Args:
        shape (float): Shape parameter of the Weibull distribution.
        scale (float, optional): Scale parameter of the Weibull distribution.
            Defaults to 1.0.

    Returns:
        Generator: A Weibull distribution sampler with given scale and shape
            parameters.

    Example:
        The following example creates a Weibull distribution sampler and 
        generates three realizations.

        ```py
        samp = weibull_sampler(1.5, 2)  # Weibull sampler with  
                                        # shape 1.5 and scale 2.0

        val1 = next(samp)   # A random value from the Weibull distribution
        val2 = next(samp)   # Another random value from the Weibull distribution
        val3 = next(samp)   # Yet another random value from the Weibull distribution
        ```
    '''
    return sampler( stats.weibull_min(shape, scale=scale) )


def beta_sampler (a, b, x0=0, x1=1):
    ''' Sampler based on a Beta distribution.

    Args:
        a (float): Shape parameter alpha.
        b (float): Shape parameter beta.
        x0 (float, optional): Lower bound of the variable definition interval.
            If omitted, it defaults to 0.
        x1 (float, optional): Upper bound of the variable definition interval.
            If omitted, it defaults to 1.

    Returns:
        Generator: A Beta distribution sampler with given shape parameters
            and definition interval boundaries.

    Example:
        The following example creates a Beta distribution sampler and 
        generates three realizations.

        ```py
        samp = beta_sampler(2, 5)  # Beta sampler with  
                                   # alpha=2, beta=5, defined in [0,1].

        val1 = next(samp)   # A random value from the Beta distribution
        val2 = next(samp)   # Another random value from the Beta distribution
        val3 = next(samp)   # Yet another random value from the Beta distribution
        ```
    '''
    return sampler( stats.beta(a, b, loc=x0, scale=x1-x0) )


def standard_gap_sampler (flange_diameter, flange_flatness_tolerance,
                          gap_angle_sampler = lognorm_sampler(100*deg, 1.0),
                          gap_shape_factor_sampler = norm_sampler(1.0, 0.15),
                          gap_height_cap = 2.0):
    ''' Sampler that generates random flange gaps according to ref. [1].

    For each realization, a gap angle is drawn from `gap_angle_sampler`, the
    gap height is then drawn from the distribution returned by
    `gap_height_distribution` for that angle, and the gap shape factor is drawn
    from `gap_shape_factor_sampler`. Realizations are discarded and redrawn
    when the gap angle is not between 0 and 180 degrees, or when the gap height
    exceeds `gap_height_cap` times the 95th percentile of its distribution.

    Args:
        flange_diameter (float): The outer diameter of the flange, in meters.
        flange_flatness_tolerance (float): The flatness tolerance in mm/mm.
        gap_angle_sampler (Generator, optional): A sampler for the gap angle,
            in radians. Defaults to a log-normal sampler with mean 100 deg and
            CoV = 1.0.
        gap_shape_factor_sampler (Generator, optional): A sampler for the gap shape factor.
            Defaults to a normal sampler with mean 1.0 and CoV = 0.15.
        gap_height_cap (float, optional): The generated gap height will be capped at
            gap_height_cap * p95, where p95 is the 95th percentile of the gap height
            distribution. If omitted, it defaults to 2.0.

    Yields:
        pyflange.flangesegments.Gap: A random gap, with random height, angle
            and shape factor.

    Example:
        The following example creates a standard gap sampler and generates
        three realizations.

        ```py
        # Gap sampler for a 7.5 m flange with flatness tolerance 1.4 mm/m
        samp = standard_gap_sampler(7.5, 0.0014)

        gap1 = next(samp)  # A random gap
        gap2 = next(samp)  # Another random gap
        gap3 = next(samp)  # Yet another random gap

        print(gap1.height, gap1.angle, gap1.shape_factor)
        ```
    '''
    from math import pi
    from .flangesegments import Gap
    from .utils import reduce_to_pi

    while True:
        gap_angle = next(gap_angle_sampler)
        gap_dist = gap_height_distribution(flange_diameter, flange_flatness_tolerance, gap_angle*flange_diameter/2)
        gap_height = gap_dist.rvs()
        if 0 < gap_angle < pi and gap_height <= gap_height_cap*gap_dist.ppf(0.95):
            yield Gap(gap_height, gap_angle, next(gap_shape_factor_sampler))



def standard_gap_sampler_ED2 (flange_diameter, flange_flatness_tolerance):
    ''' Sampler that generates random flange gaps according to ref. [3].

    For each realization, a gap angle is drawn from the distribution returned
    by `gap_angle_distribution_ED2` and the gap height is then drawn from the
    distribution returned by `gap_height_distribution_ED2` for that angle.
    The gap shape factor is always 1.0.

    Args:
        flange_diameter (float): The outer diameter of the flange, in meters.
        flange_flatness_tolerance (float): The flatness tolerance in mm/mm.

    Yields:
        pyflange.flangesegments.Gap: A random gap, with random height and angle.

    Example:
        ```py
        # ED2 gap sampler for a 7.5 m flange with flatness tolerance 1.4 mm/m
        samp = standard_gap_sampler_ED2(7.5, 0.0014)
        gap = next(samp)
        ```
    '''
    from .flangesegments import Gap

    gap_angle_dist = gap_angle_distribution_ED2(flange_diameter, flange_flatness_tolerance)
    while True:
        gap_angle = gap_angle_dist.rvs()
        gap_length = gap_angle * flange_diameter/2
        gap_height = gap_height_distribution_ED2(flange_diameter, flange_flatness_tolerance, gap_length).rvs()
        yield Gap(gap_height, gap_angle)



def product_sampler (dist1, dist2):
    ''' Sampler of the product of two independent random variables.

    Args:
        dist1 (scipy.stats.rv_frozen): Distribution of the first factor.
        dist2 (scipy.stats.rv_frozen): Distribution of the second factor.

    Yields:
        float: The product of a realization of `dist1` and a realization
            of `dist2`.

    Example:
        ```py
        from scipy.stats import norm, lognorm
        samp = product_sampler(norm(1.0, 0.1), lognorm(0.1))
        value = next(samp)
        ```
    '''
    while True:
        yield dist1.rvs() * dist2.rvs()



def clayton_copula_sampler (theta, marginal_dist1, marginal_dist2):
    ''' Sampler of two random variables joined by a Clayton copula.

    The copula realizations are generated with the conditional inverse method:
    given two independent uniform variables `u` and `w`, the second copula
    variable is

    $$v = \\left(u^{-\\theta} \\left(w^{-\\theta/(1+\\theta)} - 1\\right) + 1\\right)^{-1/\\theta}$$

    and the two variables are then mapped to the marginal distributions via
    their inverse cumulative distribution functions. The Kendall's rank
    correlation of the two variables is $\\tau = \\theta / (\\theta + 2)$.

    Args:
        theta (float): Clayton copula parameter. Must be positive.
        marginal_dist1 (scipy.stats.rv_frozen): Marginal distribution of the
            first variable.
        marginal_dist2 (scipy.stats.rv_frozen): Marginal distribution of the
            second variable.

    Yields:
        numpy.ndarray: An array `[x1, x2]` containing a joint realization of
            the two variables.

    Raises:
        ValueError: If `theta` is not positive.

    Example:
        ```py
        from scipy.stats import norm, lognorm
        samp = clayton_copula_sampler(2.0, norm(10, 1), lognorm(0.5))
        x1, x2 = next(samp)
        ```
    '''
    import numpy as np

    if not theta > 0:
        raise ValueError(f"The Clayton copula parameter must be positive, got {theta}")

    while True:
        u, w = np.random.random_sample(2)
        v = (u**(-theta) * (w**(-theta/(1+theta)) - 1) + 1)**(-1/theta)
        yield np.array([marginal_dist1.ppf(u), marginal_dist2.ppf(v)])



def clayton_gap_sampler (theta, height_marginal_dist, angle_marginal_dist, shape_factor_sampler):
    ''' Sampler that generates random flange gaps with correlated height and angle.

    The gap height and gap angle are generated as two variables joined by a
    Clayton copula (see `clayton_copula_sampler`), while the gap shape factor
    is generated independently.

    Args:
        theta (float): Clayton copula parameter. Must be positive.
        height_marginal_dist (scipy.stats.rv_frozen): Marginal distribution
            of the gap height, in meters.
        angle_marginal_dist (scipy.stats.rv_frozen): Marginal distribution
            of the gap angle, in radians.
        shape_factor_sampler (Generator): A sampler for the gap shape factor.

    Yields:
        pyflange.flangesegments.Gap: A random gap.
    '''
    from .flangesegments import Gap

    height_angle_samp = clayton_copula_sampler(theta, height_marginal_dist, angle_marginal_dist)
    while True:
        height, angle = next(height_angle_samp)
        yield Gap(height=float(height), angle=float(angle), shape_factor=next(shape_factor_sampler))



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
        ppl_sampler: Generator,         # Percentage of pretension losses sampler
        dft_sampler: Generator,         # Dry film thickness sampler
        gap_sampler: Generator,         # Gap object random sampler

        tilt_sampler: Generator = lognorm_sampler(0.1*deg, 0.50),

        E: float = 210e9,        # Young modulus of the flange
        G: float = 80.77e9,      # Shear modulus of the flange
        s_ratio: float = 1.0,    # Ratio of bottom shell thickness over s. 
                                 # Default s_botom = s.
        r: float = 0.01,         # Rounding between flange and shell
        k_shell = 'interp',      # Custom shell stiffness
        settlement_factor: float = 1.0, # Factor reducing the amount of bolt
                                        # settlement, to be taken as 1.0 for
                                        # torque tightening and 0.5 for tension
                                        # tightening.
        segment_class: type = None      # Class of the generated segments.
                                        # Default PolynomialLFlangeSegment.
    ):
    ''' A sampler that yields random pyflange PolynomialLFlangeSegment objects.

    For each realization, a random gap and a random flange tilt are drawn, and
    the bolt preload is sampled and averaged over the bolts affected by the gap
    (half of the bolts subtended by the gap angle, at least one), to account
    for the larger preload scatter at small gap angles. The returned segment
    already includes the effect of the sampled preload losses (settlement and
    plastic), so its `PPL` and `DFT` attributes are zero.

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
        ppl_sampler (Generator): Sampler of the percentage of pretension
            losses due to plastic deformation (non-dimensional, e.g. 0.03 for 3%).
        dft_sampler (Generator): Dry film thickness sampler.
        gap_sampler (Generator): Sampler of `pyflange.flangesegments.Gap`
            objects, such as `standard_gap_sampler`.
        tilt_sampler (Generator, optional): Random tilt angle sampler. Sampled
            angles are reduced to the interval [0, pi].
            Defaults to `lognorm_sampler(0.1*deg, 0.50)`.
        E (float, optional): Young modulus of the flange. Defaults to 210e9.
        G (float, optional): Shear modulus of the flange. Defaults to 80.77e9.
        s_ratio (float, optional): Ratio of bottom shell thickness over s.
            Defaults to 1.0.
        r (float, optional): Rounding between flange and shell. Defaults to 0.01.
        k_shell (str | float | None, optional): Custom shell stiffness.
            If `'interp'` (default), the stiffness is interpolated from the
            tabulated FEA results, as a function of the gap angle (see
            `pyflange.flangesegments.shell_stiffness`).
            If a number, it's used as the stiffness.
            If `None`, a simplified formula is used.
        settlement_factor (float, optional): Factor reducing the amount of bolt 
            settlement, to be taken as 1.0 for torque tightening and 0.5 for
            tension tightening. If omitted, it defaults to 1.0.
        segment_class (type, optional): The class of the generated flange
            segments: `PolynomialLFlangeSegment` or one of its subclasses,
            such as `PolynomialLFlangeSegmentED2`. If omitted, it defaults to
            `PolynomialLFlangeSegment`.

    Yields:
        pyflange.flangesegments.PolynomialLFlangeSegment: A random L-Flange segment object.

    Example:
        The following example creates a L-Flange segment sampler and generates
        three random flange segments.

        ```py
        from math import pi
        from metrum.units import mm, um, kN

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

                # bolt preload random sampler and preload losses parameters
                preload_sampler = stats.norm_sampler(2932.24*kN, 0.03), 
                ppl_sampler = stats.norm_sampler(0.03, 1/3),
                dft_sampler = stats.norm_sampler(250*um, 0.12),

                # gap object random sampler
                gap_sampler = stats.standard_gap_sampler(
                            flange_diameter = 8000*mm,
                            flange_flatness_tolerance = 0.0014,   # 1.4 mm/m
                            gap_angle_sampler = stats.lognorm_sampler(100*deg, 1.0),
                            gap_shape_factor_sampler = stats.norm_sampler(1.0, 0.15)
                        ),

                s_ratio = 1.0,              # Ratio of bottom shell thickness over s.
                settlement_factor = 0.50    # Settlement reduction for tension tightening.
            )

        fseg1 = next(fseg_samp)     # A random L-Flange segment object
        fseg2 = next(fseg_samp)     # Another random L-Flange segment object
        fseg3 = next(fseg_samp)     # Yet another random L-Flange segment object

        ```
    '''

    from math import pi
    from .flangesegments import PolynomialLFlangeSegment
    from .utils import reduce_to_pi
    import numpy as np

    segment_class = segment_class or PolynomialLFlangeSegment

    # stiffness interpolaion
    interp_shell_stiffness = _shell_stiffness_interpolator(R, s)

    # generate random flange segments
    while True:
        gap = next(gap_sampler)
        tilt = reduce_to_pi( next(tilt_sampler) )
        create_fseg = lambda Fv, PPL, DFT: segment_class(
                a=a, b=b, s=s, t=t, R=R, central_angle=central_angle, Zg=Zg, bolt=bolt,
                Do=Do, washer=washer, nut=nut, gap=gap, tilt_angle=tilt,
                E=E, G=G, s_ratio=s_ratio, r=r,
                Fv=Fv, PPL=PPL, DFT=DFT, settlement_factor=settlement_factor,
                k_shell = interp_shell_stiffness(gap.angle) if k_shell=='interp' else k_shell
            )

        # Averaging bolt preload over "50% of the bolts within the gap"
        # I am having a small sub-routine that randomly samples the preload over
        # half the gap angle and then takes the average of that for each simulation.
        # Then you get a larger scatter for smaller gap angles.
        n = max(round(gap.angle/central_angle/2), 1)  # number of bolts on which we average
        averaged_preload = np.mean([
                create_fseg(next(preload_sampler), next(ppl_sampler), next(dft_sampler)).bolt_preload
                for i in range(n)
            ])
        yield create_fseg(averaged_preload, PPL=0, DFT=0)



def _shell_stiffness_interpolator (R, s):
    # Returns a fast function gap_angle -> shell stiffness, for the given shell
    # radius and thickness, by tabulating `shell_stiffness` over 1000 angles.
    # It is also used by pyflange.io to rebuild sampled segments exactly.
    import numpy as np
    from .flangesegments import shell_stiffness

    gap_angles = np.linspace(10*deg, 180*deg, 1000)
    shell_stiffnesses = np.array([shell_stiffness(R, s, gap_angle) for gap_angle in gap_angles])
    return lambda gap_angle: np.interp(gap_angle, gap_angles, shell_stiffnesses)



def standard_markov_matrix_sampler (markov_matrix, mean_range_coeff=1.0, range_CoV=0.12):
    ''' Sampler that generates a random Markov matrix according to ref [2].

    Each realization has the same cycles, means and duration as the design
    Markov matrix, while its ranges are scaled by a single log-normally
    distributed coefficient.

    Args:
        markov_matrix (pyflange.fatigue.MarkovMatrix): The deterministic
            design Markov matrix.
        mean_range_coeff (float, optional): The mean value of the coefficient
            multiplying the load ranges of *markov_matrix*, assumed
            log-normally distributed. Defaults to 1.0.
        range_CoV (float, optional): The coefficient of variation of the
            coefficient multiplying the load ranges. Defaults to 0.12.

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



def standard_bolt_fatigue_curve_sampler (bolt_nominal_diameter, mean_stress_factor=1.0, stress_factor_CoV=0.10, DS_ref_mean=62*MPa):
    ''' Sampler that generates random bolt S-N curves, according to ref. [2].

    The reference stress range of each curve is the product of `DS_ref_mean`
    and a normally distributed random coefficient. The curves have no material
    safety factor (`gamma_M = 1.0`).

    Args:
        bolt_nominal_diameter (float): The nominal diameter of the bolt, in meters.
        mean_stress_factor (float, optional): The mean value of the fatigue class
            coefficient, assumed normally distributed. Defaults to 1.
        stress_factor_CoV (float, optional): The coefficient of variation of the
            fatigue class coefficient, assumed normally distributed. Defaults to 0.10.
        DS_ref_mean (float, optional): The mean value of the fatigue curve
            reference stress range. If unknown it could be calculated as 1.24
            times the fatigue class. If omitted, it defaults to 62 MPa.

    Yields:
        pyflange.fatigue.BoltFatigueCurve: A random bolt fatigue curve object.
    '''
    from .fatigue import BoltFatigueCurve
    stress_factor_samp = norm_sampler(mean_stress_factor, stress_factor_CoV)
    while True:
        DS_ref = DS_ref_mean * next(stress_factor_samp)
        yield BoltFatigueCurve(bolt_nominal_diameter, DS_ref, gamma_M=1.0)



def fatigue_case_sampler (fseg_samp, markov_matrix_samp, fatigue_curve_samp, allowable_damage_samp, SMF=1.0):
    ''' Sampler that generates random fatigue cases.

    Args:
        fseg_samp (Generator): A `pyflange.flangesegments.FlangeSegment` sampler,
            such as `standard_PolynomialLFlangeSegment_sampler`.
        markov_matrix_samp (Generator): A `pyflange.fatigue.MarkovMatrix` sampler.
        fatigue_curve_samp (Generator): A `pyflange.fatigue.BoltFatigueCurve` sampler.
        allowable_damage_samp (Generator): A `float` sampler that generates random allowable damages.
        SMF (float | Generator, optional): Stress Multiplication Factor. If it is a sampler, a new
            SMF will be generated for every next fatigue case. If it is a float, it will be 
            applied to every next fatigue case. If omitted, it defaults to 1.0.

    Yields:
        pyflange.fatigue.BoltFatigueAnalysis: A random `BoltFatigueAnalysis` object.
    '''
    from .fatigue import BoltFatigueAnalysis

    while True:
        yield BoltFatigueAnalysis(
            fseg = next(fseg_samp),
            flange_mkvm = next(markov_matrix_samp),
            custom_fatigue_curve = next(fatigue_curve_samp),
            allowable_damage = next(allowable_damage_samp),
            SMF = next(SMF) if isinstance(SMF, Generator) else SMF)
