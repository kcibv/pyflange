''' Fatigue calculation tools

This module defines functions and classes to support structural fatigue
calculations.

In particular, the module contains the following functions ...

- ``scf_tubular_buttweld(D, T, t, GI, GO)`` which returns the stress concentation
  factor of a butt weld between two cylindric tubulars
- ``scf_cone_buttweld(D, t, tc, alpha))`` which returns the stress concentation
  factor of a butt weld between two conic cans
- ``scf_pipe_hole(a, b)`` which returns the stress concentration factor at a hole
  in a pipe wall
- ``narrow_band_spectral_DEL(m, S, f, T, N_ref)`` which calculates the
  damage-equivalent value of a spectral load
- ``lumped_spectrum(scatter_table, w0, m)`` which returns the damage-equivalent
  spectrum from a scatter table

... and the following ``FatigueCurve`` classes:

- ``SingleSlopeFatigueCurve``
- ``DoubleSlopeFatigueCurve``

Each fatigue curve class exxposes the following methods:

- ``fatigue_curve.N(DS)`` returns the number of cycles corresponding to the
  given stress range DS
- ``fatigue_curve.DS(N)`` returns the stress range corresponding to the
  given number of cycles N
- ``fatigue_curve.damage(n, DS)`` returns the fatigue damage cumulated by
  a stress range DS repeated n times
'''

import numpy as np
from math import sqrt, pi, tan, exp, log10
from pythagoras.tools.units import mm, yr, MPa

from dataclasses import dataclass
import functools


def scf_tubular_buttweld (D, T, t, GI, GO, delta_m=4*mm):
    ''' Stress Concentration Factor for butt-welded tubulars

    This function calculates the SCF of a butt weld between two tubulars
    according to DNV-RP-C203.

    **Parameters:**

    - `D` : float
        Diameter of the connected tubulars (the two tubulars have the same diameter)

    - `T` : float
        Largest wall thickness

    - `t` : float
        Smallest wall thickness

    - `GI` : bool
        True if grinding iside is performed

    - `GO` : bool
        True if grinding outside is performed

    - `delta_m` : float
        The maximum misalignment between the two tubulars. If omitted, it
        defaults to 4 mm.

    **Return value:**
    This function returns a tuple of two values:

    - SCFInside : float
        The SCF to be applied to the inside of the weld

    - SCFOutside : float
        The SCF to be applied to the outside of the weld
    '''
    if t > T:
        t, T = T, t

    L = (T-t)*4
    delta_0 = t*0.05

    if GI == False:
        delta_total_inside = (T-t)/2 + delta_m - delta_0
    else:
        delta_total_inside = (T-t)/2 + delta_m

    if GO == False:
        delta_total_outside = (T-t)/2 - delta_m + delta_0
    else:
        delta_total_outside = (T-t)/2 - delta_m

    Beta = 1.5-1/log10(D/t) + 3/(log10(D/t)**2)
    alpha = 1.82*L/sqrt(D*t)*(1/(1+(T/t)**Beta))
    scf_inside = 1+6*(delta_total_inside/t)*(1/(1+(T/t)**Beta))*exp(-alpha)
    scf_outside = 1-6*(delta_total_outside/t)*(1/(1+(T/t)**Beta))*exp(-alpha)

    return scf_inside, scf_outside


def scf_cone_buttweld (D, t, tc, alpha):
    ''' Stress Concentration Factor for a welded cone

    scfcone returns the SCF for either the start or the bottom of a cone welded
    to a tubular, according to DNV-RP-C203.

    **Parameters:**

    - D : float
        Cylinder diameter at junction

    - t : float
        Tubular member wall thickness

    - tc : float
        Cone wall thickness

    - alpha : float
        The slope angle of the cone in radians

    **Return value:**
    This function returns a tuple of two values:

    - scf_tubular : float
        The SCF at the tubular side

    - scf_cone : float
        The SCF at the cone side.
    '''
    scf_tubular = 1+(0.6*t*sqrt(D*(t+tc))/t**2)*tan(alpha)
    scf_cone = 1+(0.6*t*sqrt(D*(t+tc))/tc**2)*tan(alpha)
    return scf_tubular, scf_cone


def scf_pipe_hole (a, b):
    ''' Stress Concentration Factor for hole in pipe wall

    This function returns the SCF occurring at a hole in the wall of a tubular.
    The hole can have an elliptical or a circular shape. The SCF is calculated
    according to Peterson's Stress Concentration Factors - Third Ed.

    **Parameters:**

    - a : float
        Horizzontal radius (half width) of the hole

    - b : float
        Vertical radius (half heght) of the hole

    **Return value:**

    - scf : float
        The stress concentration factor
    '''

    r = b**2/a
    scf = 1+2*sqrt(b/r)
    return scf


def narrow_band_spectral_DEL (S, f, m):
    ''' Narrow-band spectral damage equivalent load

    This function calculated the damage-equivalent value of a load defined by
    an energy density spectrum.

    **Parameters:**

    - `S` : n×q numpy array
        Matrix of spectral energy density values. Each column represents an
        n-dim array of spectral values corresponding to a given frequency.

    - `f` : 1×q numpy array
        Array of frequencies corresponding to the values of each column of `S`.

    - `m` : float
        Slope of the (single-sloped) fatigue curve with respect to which
        the equivalent load has to be calculated.

    **Return value:**

    - `N` : float
        Number of load cycles occurring in 1 year

    - `DEL` : 1xn numpy array
        Load range that, applied N times produces the same fatigue damage that
        the given spectrum produces in 1 year, for a single-sloped Wohler curve
        with slope m.
    '''

    # Calculate the df corresponding to each frequency
    df = np.gradient(f)

    # Calculate the spectrum moments
    m0 = S @ df            # spectrum moment of order 0
    m2 = S @ (f**2 * df)   # spectrum moment of order 2

    # Calculate the number of cycles
    Tz = (m0 / m2)**0.5         # zero-up-crossing period
    n = 1*yr / Tz               # number of cycles occurring in 1 year

    # Calculate and return the Damage Equivalent Load according to
    # Barltrop et Al. - Dynamic of fixed marine structures - 3rd ed - sec.11.3.4
    from scipy.special import gamma
    DEL = np.sqrt(8*m0) * gamma(1+m/2)**(1/m)

    # Retrun the DEL normalized to the reference number of cycles
    return n, DEL


def lumped_spectrum (scatter_table, w0, m, gamma=3.3, sigma_a=0.07, sigma_b=0.09):
    ''' Lumped spectrum around a reference frequency

    **Parameters:**

    - `scatter_table` : pandas.DataFrame
        A table of probabilities where each row corresponds to a significant
        wave height and each column to a peak period.

    - `w0` : float
        The reference angular frequency around which the spectrum is
        narrow-banded. Normally this is the first eigen frequency of the
        structure.

    - `m` : float
        The slope of the fatigue curve with respect to which the
        fatigue-equivalent sea state should be calculated.

    - *gamma* : float [optional]
        The non-dimensional peak shape parameter of the spectra defined by
        the scatter table. If omitted, it defaults to 3.3.

    - *sigma_a* : float [optional]
        The spectral width parameter for frequencies lower than 1/Tp. If
        omitted it defaults to 0.07.

    - *sigma_b* : float [optional]
        The spectral width parameter for frequencies higher than 1/Tp. If
        omitted it defaults to 0.09.

    **Return value:**

    This function returns the ``pythagoras.wtg.waves.JONSWAPSpectrum`` that produces
    alone the same fatigue damage as the entire scatter table for a Wohler
    curve with slope m, assuming that the stress spectrum is narrow-banded
    around the given angular frequency w0.
    '''
    from pythagoras.sea.waves import JONSWAPSpectrum
    Spectrum = lambda Hs, Tp: JONSWAPSpectrum(Hs, Tp, gamma, sigma_a, sigma_b)

    Hs = scatter_table.index.to_numpy()
    Tp = scatter_table.columns.to_numpy()
    P = scatter_table.to_numpy() / scatter_table.sum().sum()

    Hs_eq = np.sum(P.T @ Hs**m)**(1/m)
    S_w0 = np.array([
        [Spectrum(Hs[i], Tp[j]).Sw(w0) for j in range(len(Tp))]
        for i in range(len(Hs))
    ])
    S_w0_eq = np.sum(P * S_w0**(m/2))**(2/m)

    from scipy.optimize import fsolve
    Tp_eq = fsolve(lambda Tp_eq: abs(Spectrum(Hs_eq, Tp_eq).Sw(w0) - S_w0_eq)/S_w0_eq, np.max(Tp))[0]

    return Spectrum(Hs_eq, abs(Tp_eq))


class FatigueCurve:
    ''' A Wohler curve

    This is a base class for creating Wohler curves. It is not supposed to be
    instantiated directly.
    '''

    def N (self, DS):
        ''' Number of cycles

        Given a stress range DS, this function return the corresponding
        number of cycles that produce a fatigue failure.
        '''
        pass

    def DS (self, N):
        ''' Stress range

        Given a number of cycles, this function return the corresponding
        stress range that produce a fatigue failure.
        '''
        pass

    def damage (self, n, DS):
        ''' Fatigue damage

        Given a number of cycles n and a stress range DS, this function returns
        the dorresponding fatigue damage (D = n / N(DS)).
        '''
        return n / self.N(DS)


@dataclass
class SingleSlopeFatigueCurve (FatigueCurve):
    ''' Wohler curve with single logarithmic slope

    This class implements the FatigueCurve interface for a curve with single
    slope m.

    **Contructor parameters:**

    - `m` : float
        The logarithmic slope of the fatigue curve.

    - `DS_ref` : float
        Arbitrary reference stress range.

    - `N_ref` : float
        The number of cycles that produce failure under the stress range D_ref.

    **Attributes:**

    - `m` : float
        The logarithmic slope of the fatigue curve.

    - `DS_ref` : float
        Reference stress range.

    - `N_ref` : float
        Number of cycles at failure corresponding to the stress range D_ref.

    - `a` : float
        The Wohler curve constant (a = DS_ref**m * N_ref = DS**m * N)

    **Methods:**

    This class implements all the methods of FatigueCurve.
    '''

    m: float
    DS_ref: float
    N_ref: float


    @functools.cached_property
    def a (self):
        return self.DS_ref ** self.m * self.N_ref

    def N (self, DS):
        return self.a / DS**self.m

    def DS (self, N):
        return (self.a / N)**(1/self.m)


class MultiSlopeFatigueCurve (FatigueCurve):
    '''Multi-Slope Fatigue Curve

    This class is a FatigueCurve with multiple slopes.
    It takes any number of SingleSlopeFatigueCurve objects as arguments.
    '''

    def __init__ (self, *fatigue_curves):
        self.curves = fatigue_curves

    def N (self, DS):
        return max([curve.N(DS) for curve in self.curves])

    def DS (self, N):
        return max([curve.DS(N) for curve in self.curves])


class DoubleSlopeFatigueCurve (MultiSlopeFatigueCurve):
    ''' Wohler curve with double logarithmic slope

    This class implements the FatigueCurve interface for a curve with two
    slopes m1 and m2.

    **Contructor parameters:**

    - `m1` : float
        The logarithmic slope of the lower cycle values.

    - `m2` : float
        The logarithmic slope of the higher cycle values.

    - `DS12` : float
        The stress range where the two branches of the curve meet.

    - `N12` : float
        The number of cycles to failure corresponding to DS12.

    **Attributes:**

    - `m1` : float
        The logarithmic slope of the lower cycle values.

    - `m2` : float
        The logarithmic slope of the higher cycle values.

    - `DS12` : float
        The stress range where the two branches of the curve meet.

    - `N12` : float
        The number of cycles to failure corresponding to DS12.

    **Methods:**

    This class implements all the methods of FatigueCurve.
    '''

    def __init__ (self, m1, m2, DS12, N12):
        curve1 = SingleSlopeFatigueCurve(m1, DS12, N12)
        curve2 = SingleSlopeFatigueCurve(m2, DS12, N12)
        super().__init__(curve1, curve2)


dnv_fatigue_curves = {

    # S-N curves in air
    "B1 (A)": DoubleSlopeFatigueCurve(4, 5, 106.97*MPa, 1e7),
    "B2 (A)": DoubleSlopeFatigueCurve(4, 5,  93.59*MPa, 1e7),
    "C (A)" : DoubleSlopeFatigueCurve(3, 5,  73.10*MPa, 1e7),
    "C1 (A)": DoubleSlopeFatigueCurve(3, 5,  65.50*MPa, 1e7),
    "C2 (A)": DoubleSlopeFatigueCurve(3, 5,  58.48*MPa, 1e7),
    "D (A)" : DoubleSlopeFatigueCurve(3, 5,  52.63*MPa, 1e7),
    "E (A)" : DoubleSlopeFatigueCurve(3, 5,  46.78*MPa, 1e7),
    "F (A)" : DoubleSlopeFatigueCurve(3, 5,  41.52*MPa, 1e7),
    "F1 (A)": DoubleSlopeFatigueCurve(3, 5,  36.84*MPa, 1e7),
    "F3 (A)": DoubleSlopeFatigueCurve(3, 5,  32.75*MPa, 1e7),
    "G (A)" : DoubleSlopeFatigueCurve(3, 5,  29.24*MPa, 1e7),
    "W1 (A)": DoubleSlopeFatigueCurve(3, 5,  26.32*MPa, 1e7),
    "W2 (A)": DoubleSlopeFatigueCurve(3, 5,  23.39*MPa, 1e7),
    "W3 (A)": DoubleSlopeFatigueCurve(3, 5,  21.05*MPa, 1e7),

    # S-N curves in seawater with cathodic protection
    "B1 (W)": DoubleSlopeFatigueCurve(4, 5, 106.97*MPa*10**(0.2), 1e6),
    "B2 (W)": DoubleSlopeFatigueCurve(4, 5,  93.59*MPa*10**(0.2), 1e6),
    "C (W)" : DoubleSlopeFatigueCurve(3, 5,  73.10*MPa*10**(0.2), 1e6),
    "C1 (W)": DoubleSlopeFatigueCurve(3, 5,  65.50*MPa*10**(0.2), 1e6),
    "C2 (W)": DoubleSlopeFatigueCurve(3, 5,  58.48*MPa*10**(0.2), 1e6),
    "D (W)" : DoubleSlopeFatigueCurve(3, 5,  52.63*MPa*10**(0.2), 1e6),
    "E (W)" : DoubleSlopeFatigueCurve(3, 5,  46.78*MPa*10**(0.2), 1e6),
    "F (W)" : DoubleSlopeFatigueCurve(3, 5,  41.52*MPa*10**(0.2), 1e6),
    "F1 (W)": DoubleSlopeFatigueCurve(3, 5,  36.84*MPa*10**(0.2), 1e6),
    "F3 (W)": DoubleSlopeFatigueCurve(3, 5,  32.75*MPa*10**(0.2), 1e6),
    "G (W)" : DoubleSlopeFatigueCurve(3, 5,  29.24*MPa*10**(0.2), 1e6),
    "W1 (W)": DoubleSlopeFatigueCurve(3, 5,  26.32*MPa*10**(0.2), 1e6),
    "W2 (W)": DoubleSlopeFatigueCurve(3, 5,  23.39*MPa*10**(0.2), 1e6),
    "W3 (W)": DoubleSlopeFatigueCurve(3, 5,  21.05*MPa*10**(0.2), 1e6),

    # S-N curves for free corrosion in plated structures
    "B1 (F)": SingleSlopeFatigueCurve(3, 10**((12.436-7)/3)*MPa, 1e7),
    "B2 (F)": SingleSlopeFatigueCurve(3, 10**((12.262-7)/3)*MPa, 1e7),
    "C (F)" : SingleSlopeFatigueCurve(3, 10**((12.115-7)/3)*MPa, 1e7),
    "C1 (F)": SingleSlopeFatigueCurve(3, 10**((11.972-7)/3)*MPa, 1e7),
    "C2 (F)": SingleSlopeFatigueCurve(3, 10**((11.824-7)/3)*MPa, 1e7),
    "D (F)" : SingleSlopeFatigueCurve(3, 10**((11.687-7)/3)*MPa, 1e7),
    "E (F)" : SingleSlopeFatigueCurve(3, 10**((11.533-7)/3)*MPa, 1e7),
    "F (F)" : SingleSlopeFatigueCurve(3, 10**((11.378-7)/3)*MPa, 1e7),
    "F1 (F)": SingleSlopeFatigueCurve(3, 10**((11.222-7)/3)*MPa, 1e7),
    "F3 (F)": SingleSlopeFatigueCurve(3, 10**((11.068-7)/3)*MPa, 1e7),
    "G (F)" : SingleSlopeFatigueCurve(3, 10**((10.921-7)/3)*MPa, 1e7),
    "W1 (F)": SingleSlopeFatigueCurve(3, 10**((10.784-7)/3)*MPa, 1e7),
    "W2 (F)": SingleSlopeFatigueCurve(3, 10**((10.630-7)/3)*MPa, 1e7),
    "W3 (F)": SingleSlopeFatigueCurve(3, 10**((10.493-7)/3)*MPa, 1e7)
}
