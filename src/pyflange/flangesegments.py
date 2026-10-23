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

'''Mechanical models of single-bolt flange segments.

This module contains `FlangeSegment` classes, which model the mechanical
behavior of a flange sector containing one bolt only. The module provides:

- `Gap`: description of the imperfection gap between two flange faces.
- `FlangeSegment`: abstract base class defining the transfer-function interface
  (bolt axial force and bolt bending moment as functions of the shell pull force).
- `PolynomialFlangeSegment`: abstract base class implementing the polynomial
  transfer functions of refs. [1] and [2].
- `PolynomialLFlangeSegment` and `PolynomialLFlangeSegmentED2`: polynomial
  L-flange segment models (the latter implementing the modifications expected
  in IEC 61400-6 ED2).
- `PolynomialTFlangeSegment` and `PolynomialTFlangeSegmentED2`: polynomial
  T-flange segment models (the latter implementing the modifications expected
  in IEC 61400-6 ED2).
- `bolt_markov_matrix`: conversion of a flange bending-moment Markov matrix into
  a bolt stress Markov matrix, via the segment transfer functions.
- `shell_stiffness`: shell stiffness interpolated from tabulated FEA results.

The only types of FlangeSegment currently available are an L-Flange segment
and a T-Flange segment, implementing a polynomial relation between tower shell
force and bolt force / bolt moment. Nonetheless, this module has been structured
to be easily extensible with other types of FlangeSegment model, such as
Multilinear (Petersen) L-Flanges, Multilinear T-Flanges, etc.

Units:
    All parameters must be supplied in a consistent system of units, e.g. meters
    for lengths, newtons for forces and pascals for stresses. Angles are
    expressed in radians.

Example:
    The following example creates an L-flange segment and evaluates the bolt
    axial force for a given shell pull force.

    ```python
    from math import pi
    from pyflange.bolts import StandardMetricBolt, RoundNut
    from pyflange.flangesegments import Gap, PolynomialLFlangeSegment

    fseg = PolynomialLFlangeSegment(
        a=0.150, b=0.122, s=0.054, t=0.172, R=4.0,
        central_angle=2*pi/156,
        Zg=-115e3,                              # dead weight per bolt (compression)
        bolt=StandardMetricBolt("M80", "10.9", shank_length=0.160, stud=True),
        Fv=2.9e6,                               # bolt preload
        Do=0.086, washer=None, nut=RoundNut("M80"),
        gap=Gap(height=0.002, angle=pi/6))

    Fs = fseg.bolt_axial_force(500e3)           # bolt force for a shell pull of 500 kN
    ```

References:
    The reference numbers are the same in all the PyFlange modules (see the
    [references page](../references.md)). This module cites:

    - `[1]` Seidel, M.: IEC 61400-6 AMD1 Background document. Fatigue design
      for ring flange connections in wind turbine support structures. Rev00,
      Siemens Gamesa, 2024.
    - `[2]` IEC 61400-6:2020+AMD1:2025 CSV. Wind energy generation systems -
      Part 6: Tower and foundation design requirements. Consolidated version,
      2025. Equation numbers G.xx refer to its Annex G.
    - `[3]` Petersen, C.: Nachweis der Betriebsfestigkeit exzentrisch
      beanspruchter Ringflanschverbindungen. Stahlbau 67 (1998), pp. 191-203.
      <https://doi.org/10.1002/stab.199800690>
    - `[4]` Petersen, C.: Stahlbau. 4th edition, Springer Vieweg, 2012.
      <https://doi.org/10.1007/978-3-8348-8610-1>
    - `[5]` Seidel, M.: IEC 61400-6 AMD1 & ED2 Background document. Fatigue
      design for ring flange connections in wind turbine support structures.
      Rev01, draft 2026-08-04, Siemens Gamesa. The changes of the polynomial
      model from AMD1 to ED2 are summarised in its section 9.7.
    - `[6]` Tobinaga, I.; Ishihara, T.: A study of action point correction
      factor for L-type flanges of wind turbine towers. Wind Energy 21 (2018),
      pp. 801-806. <https://doi.org/10.1002/we.2193>
'''



from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cache, cached_property

from .utils import Logger, log_data
logger = Logger(__name__)

import numpy as np

from .bolts import Bolt, Washer, Nut

from math import pi

import pandas as pd
#from.fatigue import SNCurve



@dataclass
class Gap:
    '''A flange gap due to imperfections.

    The design gap height of IEC 61400-6 AMD1 is the 95% quantile of the
    log-normal distribution of ref. [2], section 6.7.5.2 (see
    `pyflange.stats.gap_height_distribution`). The gap angles to be checked
    are 30°, 60°, 90° and 120° (ref. [2], section 6.7.5.3).

    Attributes:
        height (float): Maximum longitudinal gap height.
        angle (float): Angle subtended by the gap arc from the flange center,
            in radians.
        shape_factor (float): Factor accounting for a shape different than
            sinusoidal. Defaults to 1.0 (sinusoidal gap).
    '''

    height       : float         # maximum longitudinal gap height
    angle        : float         # angle subtended by the gap arc from the flange center
    shape_factor : float = 1.0   # Factor accounting for a shape different than sinusoidal



class FlangeSegment (ABC):
    '''Abstract base class for a flange segment.

    This class is meant to be extended and not to be instantiated directly.
    Each FlangeSegment child class must implement the two transfer-function
    methods: `.bolt_axial_force(Z)` and `.bolt_bending_moment(Z)`.
    '''

    @abstractmethod
    def bolt_axial_force (self, shell_pull):
        '''Calculates the bolt axial force generated by a shell pull force.

        Args:
            shell_pull (float or numpy.ndarray): The shell pull force Z.

        Returns:
            float or numpy.ndarray: The resulting bolt axial force.
        '''
        pass

    @abstractmethod
    def bolt_bending_moment (self, shell_pull):
        '''Calculates the bolt bending moment generated by a shell pull force.

        Args:
            shell_pull (float or numpy.ndarray): The shell pull force Z.

        Returns:
            float or numpy.ndarray: The resulting bolt bending moment.
        '''
        pass



class PolynomialFlangeSegment (FlangeSegment):
    '''A generic FlangeSegment with polynomial transfer functions.

    This class implements a polynomial relation between tower shell force Z and
    bolt axial force Fs or bolt bending moment Ms. It is not meant to be
    instantiated directly but to be subclassed.

    The polynomial transfer functions are defined by 4 reference points:

    - P1: Flange segment at rest.
    - P2: Flange segment at ultimate tensile limit state.
    - P3: Flange segment in small tensile deformation.
    - P4: Gap closure state.

    Child classes must provide the coordinates of these points (shell force,
    bolt force and bolt moment) through the abstract properties of this class.
    The transfer functions are a quadratic polynomial through P1, P2 and P3 in
    the tensile domain (above the shell force at rest), a quadratic polynomial
    between P4 and P1 in the compressive domain, and a constant below P4.
    Child classes may override the tensile polynomials when the reference
    defines them differently (e.g. the cubic bolt moment of T-flanges).
    '''

    def bolt_axial_force (self, shell_pull):
        '''Calculates bolt axial force from tower shell force using a polynomial.

        The relation is defined as a polynomial function based on
        ref. [2], section G.2.2, eq. G.32.

        Args:
            shell_pull (float or numpy.ndarray): The tower shell force Z.

        Returns:
            float or numpy.ndarray: The corresponding bolt axial force Fs.
        '''

        # Retrueve the polynomial branch in the tensile domain and the
        # polynomial branch in the compressive domain. The evaluation
        # of both polynomials is delegated to a separate getter for
        # the sake of separation of concerns and code readability.
        tens_poli = self._tensile_force_polynomial
        comp_poli = self._compressive_force_polynomial

        # The tensile polynomial is defined between Z1 and Z2, while
        # the compressive prolynomila between Zmin and Z1
        Z1, Z2 = tens_poli.domain
        Zmin = comp_poli.domain[0]

        # Calculate a vector with the same size as shell_pull, having the Fs value
        # in tensile domain where shell_value is greater than Z1 and 0 in all
        # the other cases.
        Fs_tens = tens_poli(shell_pull) * (shell_pull > Z1)

        # Calculate a vector with the same size as shell_pull, having the Fs value
        # in compressive domain where shell_value is between Zmin and Z1 and 0 in all
        # the other cases.
        Fs_comp = comp_poli(shell_pull) * (shell_pull > Zmin) * (shell_pull <= Z1)

        # Calculate a vector with the same size as shell_pull, having the Fs(Zmin)
        # value where shell_value is lower than or equal to zmin and 0 in all
        # the other cases.
        Fs_min  = comp_poli(Zmin) * (shell_pull <= Zmin)

        # Compose the three branches together and return them.
        return Fs_tens + Fs_comp + Fs_min


    def bolt_bending_moment (self, shell_pull):
        '''Calculates bolt bending moment from tower shell force using a polynomial.

        The relation is defined as a polynomial function based on
        ref. [2], section G.2.3, eqs. G.59 (L-flanges) and G.61 (T-flanges).

        Args:
            shell_pull (float or numpy.ndarray): The tower shell force Z.

        Returns:
            float or numpy.ndarray: The corresponding bolt bending moment Ms.
        '''

        # Retrueve the polynomial branch in the tensile domain and the
        # polynomial branch in the compressive domain. The evaluation
        # of both polynomials is delegated to a separate getter for
        # the sake of separation of concerns and code readability.
        tens_poli = self._tensile_moment_polynomial
        comp_poli = self._compressive_moment_polynomial

        # The tensile polynomial is defined between Z1 and Z2, while
        # the compressive prolynomila between Zmin and Z1
        Z1, Z2 = tens_poli.domain
        Zmin = comp_poli.domain[0]

        # Calculate a vector with the same size as shell_pull, having the Ms values
        # in tensile domain where shell_value is greater than Z1 and 0 in all
        # the other cases.
        Ms_tens = tens_poli(shell_pull) * (shell_pull > Z1)

        # Calculate a vector with the same size as shell_pull, having the Ms values
        # in compressive domain where shell_value is between Zmin and Z1 and 0 in all
        # the other cases.
        Ms_comp = comp_poli(shell_pull) * (shell_pull > Zmin) * (shell_pull <= Z1)

        # Calculate a vector with the same size as shell_pull, having the Ms(Zmin)
        # value where shell_value is lower than or equal to zmin and 0 in all
        # the other cases.
        Ms_min  = comp_poli(Zmin) * (shell_pull <= Zmin)

        # Compose the three branches together and return them.
        return Ms_tens + Ms_comp + Ms_min


    @cached_property
    def _tensile_force_polynomial (self):
        '''Polynomial Fs(Z) for the tensile domain.

        This getter returns the polynomial Fs(Z) defined by refs. [1] and [2]
        between point 1 (flange segment at rest) and point 2 (ultimate
        tensile limit state of the bolt).

        Returns:
            numpy.polynomial.polynomial.Polynomial: The polynomial for the tensile domain.
        '''

        from numpy.polynomial.polynomial import Polynomial
        from numpy.linalg import inv

        # The polynomial function Fs(Z) passes through the following points
        # The point getters are implementation-specific and they should be
        # defined by the child class (L-Flange class or T-Flange class).
        Z1 = self.shell_force_at_rest
        Fs1 = self.bolt_force_at_rest

        Z2 = self.shell_force_at_tensile_ULS
        Fs2 = self.bolt_force_at_tensile_ULS

        Z3 = self.shell_force_at_small_displacement
        Fs3 = self.bolt_force_at_small_displacement

        # The three fitting conditions Fs(Z1)=Fs1, Fs(Z2)=Fs2 and Fs(Z3)=Fs3
        # are written in matrix form as: Z*A=F whenre A is the coefficients
        # vector and F and Z are defined below:
        F = np.array([Fs1, Fs2, Fs3])
        Z = np.array([
            [1, Z1, Z1**2],
            [1, Z2, Z2**2],
            [1, Z3, Z3**2]])

        # The coefficients are therefore:
        a0, a1, a2 = inv(Z) @ F

        # Create and return the polynomial for the tensile domain [Z1, Z2]
        return Polynomial(
            (a0, a1, a2),       # coefficients
            domain=(Z1, Z2),    # the polynomial is defined between Z1 and Z2
            window=(Z1, Z2),
            symbol="Z")         # used to render the independent variable in serializations


    @cached_property
    def _compressive_force_polynomial (self):
        '''Polynomial Fs(Z) for the compressive domain.

        This getter returns the polynomial Fs(Z) defined by refs. [1] and [2]
        between point 1 (flange segment at rest) and point 4 (gap closed).

        Returns:
            numpy.polynomial.polynomial.Polynomial: The polynomial for the compressive domain.
        '''

        from numpy.polynomial.polynomial import Polynomial

        # Retrieve origin of the compressive polynomial
        Z1 = self.shell_force_at_rest
        Fs1 = self.bolt_force_at_rest

        # Retrieve slope of the compressive polynomial at point (Z1, Fs1)
        dtp = self._tensile_force_polynomial.deriv()    # first derivative of the tensile polynomial
        X1 = dtp(Z1)                                    # slope of the tensile polinomial at Z=Z1

        # Value of Z below which the compressive polynomial
        # becomes practically constant
        Zmin = self.shell_force_at_closed_gap

        # Compressive polynomial coefficients
        c2 = 0.5 * X1 / (Z1 - Zmin)
        c1 = X1 - 2*c2*Z1
        c0 = Fs1 - c1*Z1 - c2*Z1**2

        # Create and return the polynomial for the tensile domain [Zmin, Z1]
        return Polynomial(
            (c0, c1, c2),       # coefficients
            domain=(Zmin, Z1),  # the polynomial is defined between Zmin and Z1
            window=(Zmin, Z1),
            symbol="Z")         # used to render the independent variable in serializations


    @cached_property
    def _tensile_moment_polynomial (self):
        '''Polynomial Ms(Z) for the tensile domain.

        This getter returns the polynomial Ms(Z) defined by refs. [1] and [2]
        between point 1 (flange segment at rest) and point 2 (ultimate
        tensile limit state of the bolt).

        Returns:
            numpy.polynomial.polynomial.Polynomial: The polynomial for the tensile domain.
        '''

        from numpy.polynomial.polynomial import Polynomial
        from numpy.linalg import inv

        # The polynomial function Fs(Z) passes through the following points
        # The point getters are implementation-specific and they should be
        # defined by the child class (L-Flange class or T-Flange class).
        Z1 = self.shell_force_at_rest
        Ms1 = self.bolt_moment_at_rest

        Z2 = self.shell_force_at_tensile_ULS
        Ms2 = self.bolt_moment_at_tensile_ULS

        Z3 = self.shell_force_at_small_displacement
        Ms3 = self.bolt_moment_at_small_displacement

        # The three fitting conditions Fs(Z1)=Fs1, Fs(Z2)=Fs2 and Fs(Z3)=Fs3
        # are written in matrix form as: Z*A=F whenre A is the coefficients
        # vector and F and Z are defined below:
        M = np.array([Ms1, Ms2, Ms3])
        Z = np.array([
            [1, Z1, Z1**2],
            [1, Z2, Z2**2],
            [1, Z3, Z3**2]])

        # The coefficients are therefore:
        a0, a1, a2 = inv(Z) @ M

        # Create and return the polynomial for the tensile domain [Z1, Z2]
        return Polynomial(
            (a0, a1, a2),       # coefficients
            domain=(Z1, Z2),    # the polynomial is defined between Z1 and Z2
            window=(Z1, Z2),
            symbol="Z")         # used to render the independent variable in serializations


    @cached_property
    def _compressive_moment_polynomial (self):
        '''Polynomial Ms(Z) for the compressive domain.

        This getter returns the polynomial Ms(Z) defined by refs. [1] and [2]
        between point 1 (flange segment at rest) and point 4 (gap closed).

        Returns:
            numpy.polynomial.polynomial.Polynomial: The polynomial for the compressive domain.
        '''

        from numpy.polynomial.polynomial import Polynomial

        # Retrieve origin of the compressive polynomial
        Z1 = self.shell_force_at_rest
        Ms1 = self.bolt_moment_at_rest

        # Retrieve slope of the compressive polynomial at point (Z1, Fs1)
        dtp = self._tensile_moment_polynomial.deriv()    # first derivative of the tensile polynomial
        X1 = dtp(Z1)                                    # slope of the tensile polinomial at Z=Z1

        # Value of Z below which the compressive polynomial
        # becomes practically constant
        Zmin = self.shell_force_at_closed_gap

        # Compressive polynomial coefficients
        c2 = 0.5 * X1 / (Z1 - Zmin)
        c1 = X1 - 2*c2*Z1
        c0 = Ms1 - c1*Z1 - c2*Z1**2

        # Create and return the polynomial for the tensile domain [Zmin, Z1]
        return Polynomial(
            (c0, c1, c2),       # coefficients
            domain=(Zmin, Z1),  # the polynomial is defined between Zmin and Z1
            window=(Zmin, Z1),
            symbol="Z")         # used to render the independent variable in serializations


    @property
    @abstractmethod
    def shell_force_at_rest (self):
        '''Shell force when no external loads are applied (Z1).'''
        pass


    @property
    @abstractmethod
    def bolt_force_at_rest (self):
        '''Bolt axial force when no external loads are applied (Fs1).'''
        pass


    @property
    @abstractmethod
    def bolt_moment_at_rest (self):
        '''Bolt bending moment when no external loads are applied (Ms1).'''
        pass


    @property
    @abstractmethod
    def shell_force_at_small_displacement (self):
        '''Intermediate tower shell force (Z3).'''
        pass


    @property
    @abstractmethod
    def bolt_force_at_small_displacement (self):
        '''Intermediate bolt axial force (Fs3).'''
        pass


    @property
    @abstractmethod
    def bolt_moment_at_small_displacement (self):
        '''Intermediate bolt bending moment (Ms3).'''
        pass


    @property
    @abstractmethod
    def shell_force_at_tensile_ULS (self):
        '''Ultimate Limit State tower shell force (Z2).'''
        pass


    @property
    @abstractmethod
    def bolt_force_at_tensile_ULS (self):
        '''Bolt axial force at tensile failure (Fs2).'''
        pass


    @property
    @abstractmethod
    def bolt_moment_at_tensile_ULS (self):
        '''Bolt bending moment at tensile failure (Ms2).'''
        pass


    @property
    @abstractmethod
    def shell_force_at_closed_gap (self):
        '''Shell force necessary to completely close the imperfection gap (Z4).'''
        pass



@dataclass
class PolynomialLFlangeSegment (PolynomialFlangeSegment):
    '''An L-Flange implementation of a `PolynomialFlangeSegment`.

    This class models an L-flange segment using the polynomial model from
    ref. [1] and [2]. It defines the reference points for the polynomial
    functions inherited from the parent class.

    The transfer functions are available through the `bolt_axial_force` and
    `bolt_bending_moment` methods, which accept either a scalar or a NumPy
    array of shell pull forces.

    Attributes:
        a (float): Distance from inner flange face to bolt hole center.
        b (float): Distance from bolt hole center to shell center-line.
        s (float): Shell thickness.
        t (float): Flange thickness.
        R (float): Shell outer curvature radius.
        central_angle (float): Angle subtended by the flange segment arc.
        Zg (float): Load on the shell at rest (e.g., dead weight). Negative for compression.
        bolt (Bolt): The bolt in the flange segment.
        Fv (float): Design preload value before losses.
        Do (float): Bolt hole diameter.
        washer (Washer): The washer used. None if no washer.
        nut (Nut): The nut used.
        gap (Gap): The flange gap.
        tilt_angle (float): Flange tilt angle in radians.
        E (float): Young's modulus of the flange material.
        G (float): Shear modulus of the flange material.
        s_ratio (float): Ratio of bottom shell thickness to `s`.
        r (float): Radius of the rounding between shell and flange.
        k_shell (float, callable, None): Shell stiffness [N/m/m]. Can be a value, a
            function taking the flange segment as argument and returning the
            stiffness, or `None` to use the simplified formula of ref. [2], eq. G.43.
        PPL (float): Percentage of bolt pretension loss due to plastic strain,
            expressed as a fraction of `Fv` (e.g. 0.03 for 3%). Defaults to 0.0.
        DFT (float): Maximum allowable dry film thickness of coating applied on the contact
            area between washers and flange [m]. Defaults to 0.0 m.
        settlement_factor (float): Factor reducing the amount of bolt settlement, to be
            taken as 1.0 for torque tightening and 0.5 for tension tightening. Defaults to 1.0.
    '''

    a: float        # distance between inner face of the flange and center of the bolt hole
    b: float        # distance between center of the bolt hole and center-line of the shell
    s: float        # shell thickness
    t: float        # flange thickness
    R: float        # shell outer curvature radius
    central_angle: float     # angle subtended by the flange segment

    Zg: float       # load applied to the flange segment shell at rest (normally dead weight
                    # of tower + RNA, divided by the number of bolts). Negative if compression.

    bolt: Bolt      # Bolt object representing the flange segment bolt
    Fv: float       # applied bolt preload before pretension losses

    Do: float       # Bolt hole diameter
    washer: Washer  # Bolt washer
    nut: Nut        # Bolt nut

    gap: Gap        # Object representing the flange gap

    tilt_angle: float = 0.0         # angle of flange tilt

    E: float = 210e9        # Young modulus of the flange
    G: float = 80.77e9      # Shear modulus of the flange
    s_ratio: float = 1.0    # Ratio of bottom shell thickness over s. Default s_botom = s.
    r: float = 0.01         # Rounding between flange and shell
    
    k_shell: float = None       # optional initial shell stiffness in [N/m/m].

    PPL: float = 0.0      # Percentage of pretension losses due to plastic strain

    DFT: float = 0.0        # Maximum allowable dry film thickness of coatings applied in the 
                            # contact area between washers and flange [m]

    settlement_factor: float = 1.0   # Factor reducing the amount of bolt settlement, to be
                                     # taken as 1.0 for torque tightening and 0.5 for tension tightening.


    def failure_mode (self, fy_sh, fy_fl, gamma_0=1.1):
        '''Determines the failure mode of the flange.

        The ultimate shell force is evaluated for each of the failure modes
        A, B, D and E of ref. [2], and the governing mode is the one with the
        lowest ultimate shell force.

        Args:
            fy_sh (float): Yield stress of the shell material.
            fy_fl (float): Yield stress of the flange material.
            gamma_0 (float, optional): Material factor. Defaults to 1.1.

        Returns:
            tuple: A tuple containing:

                - str: The failure mode ("A", "B", "D", or "E").
                - list: A list of the shell force limits for each failure mode [Zu_A, Zu_B, Zu_D, Zu_E].
        '''

        from scipy.optimize import fsolve

        fd_sh = fy_sh / gamma_0
        fd_fl = fy_fl / gamma_0
        
        Dw = self.washer.outer_diameter if self.washer else self.nut.bearing_diameter
        b_red = self.b-self.s/2-0.8*self.r

        # Segment width at the middle surface of the shell (ref. [2], eq. G.34)
        # and segment width reduced by the bolt hole (ref. [2], eq. G.7)
        c = self.central_angle * (self.R - self.s/2)
        c_hole = c - self.Do

        # Failure mode A
        F_tRd = self.bolt.ultimate_tensile_capacity()   # Bolt ultimate tensile capacity
        Zu_A = F_tRd                                    # Ultimate tower shell force for failure mode A

        # Shell cross-section ultimate capacities
        Nu_sh = fd_sh * c * self.s                 # Pure axial ultimate capacity
        Mu_sh = fd_sh * c * self.s**2 / 4          # Pure bending ultimate capacity
        MNu_sh = lambda N: Mu_sh * (1 - (N / Nu_sh)**2) if N < Nu_sh else 0 # Bending ultimate capacity, concurrent with axial force

        # Flange cross-section ultimate capacity, ref. [2], eqs. G.8, G.10 and G.11
        Vu_fl = fd_fl * c * self.t / 3**0.5                # Pure shear ultimate capacity
        Mu_fl = fd_fl * c * self.t**2 / 4                  # Pure bending ultimate capacity
        MVu_fl = lambda V: Mu_fl * (1 - (V / Vu_fl)**2)**0.5 if V < Vu_fl else 0    # Bending ultimate capacity, concurrent with shear

        # Iterative calculation of the ultimate shell force, given the two functions:
        # - Zu_sh(Mu_pl3), returning the shell force at shell axial ULS, for a given value of Mu_pl3
        # - Zu_fl(Mu_pl3), returning the shell force at slange shear ULS, for a given value of Mu_pl3
        # Ref. [2], eq. G.11.
        def find_Zu (Zu_sh, Zu_fl):
            target_Zu_sh = fsolve(lambda Z: Zu_sh(MNu_sh(Z)) - Z, Zu_sh(Mu_sh))[0]
            target_Zu_fl = fsolve(lambda Z: Zu_fl(MVu_fl(Z)) - Z, Zu_fl(Mu_fl))[0]
            return min(target_Zu_sh, target_Zu_fl)

        # Failure mode B, ref. [2], eqs. G.3 and G.12
        Zu_sh_B = lambda Mu_pl3: (F_tRd * self._a_red + Mu_pl3) / (self._a_red + self.b)
        Zu_fl_B = lambda Mu_pl3: (F_tRd * self._a_red + Mu_pl3) / (self._a_red + b_red)
        Zu_B = find_Zu(Zu_sh_B, Zu_fl_B)

        # Failure mode D, ref. [2], eqs. G.4, G.7 and G.9
        Mu_pl2 = fd_fl * c_hole * self.t**2 / 4
        DMu_pl2 = F_tRd/2 * (self.Do/2 + Dw/2)/2
        Zu_sh_D = lambda Mu_pl3: (Mu_pl2 + DMu_pl2 + Mu_pl3) / self.b
        Zu_fl_D = lambda Mu_pl3: (Mu_pl2 + DMu_pl2 + Mu_pl3) / b_red
        Zu_D = find_Zu(Zu_sh_D, Zu_fl_D) 

        # Failure mode E, ref. [2], eq. G.5
        b_E = self.b - (self.Do+Dw)/4
        b_E_red = b_E - self.s/2 - 0.8*self.r
        Zu_sh_E = lambda Mu_pl3: (Mu_fl + Mu_pl3) / b_E
        Zu_fl_E = lambda Mu_pl3: (Mu_fl + Mu_pl3) / b_E_red
        Zu_E = find_Zu(Zu_sh_E, Zu_fl_E) 

        # Determine the governing failure mode
        Zu_min = min(Zu_A, Zu_B, Zu_D, Zu_E)

        if Zu_min == Zu_A:
            return "A",[Zu_A, Zu_B, Zu_D, Zu_E]
        elif Zu_min == Zu_B:
            return "B",[Zu_A, Zu_B, Zu_D, Zu_E]
        elif Zu_min == Zu_D:
            return "D",[Zu_A, Zu_B, Zu_D, Zu_E]
        elif Zu_min == Zu_E:
            return "E",[Zu_A, Zu_B, Zu_D, Zu_E]


    def validate (self, fy_sh, fy_fl, gamma_0=1.1):
        '''Verifies the applicability of the polynomial model.

        Checks if this L-Flange Segment fails according to failure mode B.

        Args:
            fy_sh (float): Yield stress of the shell material.
            fy_fl (float): Yield stress of the flange material.
            gamma_0 (float, optional): Material factor. Defaults to 1.1.

        Raises:
            ValueError: If the failure mode is not "B".
        '''

        failure_mode, failure_shellforces = self.failure_mode(fy_sh, fy_fl, gamma_0=gamma_0)

        if failure_mode != "B":
            raise ValueError(f"Failure mode {failure_mode} detected, while only failure mode B is supported.")


    @cached_property
    def _bolt_axial_stiffness (self):
        '''The axial stiffness of the bolt.'''
        tw = self.washer.thickness if self.washer else 0
        return self.bolt.axial_stiffness(2*self.t + 2*tw)


    @cached_property
    def _bolt_bending_stiffness (self):
        '''The bending stiffness of the bolt.'''
        tw = self.washer.thickness if self.washer else 0
        return self.bolt.bending_stiffness(2*self.t + 2*tw)


    @cached_property
    def _flange_axial_stiffness (self):
        '''The axial stiffness of the flange in the thickness direction.'''
        # Stiffnes of flange w.r.t. compression in thickness direction,
        # when no gap is present. Calculated according to ref. [3] and [4].
        from math import pi, inf
        Dw = self.washer.outer_diameter if self.washer else self.nut.bearing_diameter
        Do = self.Do
        h = self.t * 2
        A = pi * ((Dw + h/10)**2 - Do**2) / 4
        Kf = self.E * A / h
        Kw = self.washer.axial_stiffness if self.washer else inf
        return 1 / (1/Kf + 2 * 1/Kw)


    def _bolt_moment (self, Z, Fs):
        '''Calculates the bolt bending moment, according to ref. [2], eq. G.59.

        Args:
            Z (float): Shell pull force.
            Fs (float): Bolt axial force.

        Returns:
            float: The bolt bending moment.
        '''
        s_avg = self.s * (1 + self.s_ratio) / 2
        c = self.central_angle * (self.R - s_avg/2)
        I_tg = c * self.t**3 / 12
        ak = self._stiffness_correction_factor
        bolt_rotation = Z*self.b*self._a_star / (3*self.E*I_tg*ak) + (Fs - self.bolt_preload) / (2*self._a_star*self._bolt_axial_stiffness)

        log_data(self, I_tg=I_tg)

        return bolt_rotation * 2*self._bolt_bending_stiffness


    @cached_property
    def bolt_preload (self):
        '''float: Bolt preload after settlement and plastic losses.

        The preload `Fv` is reduced by the settlement loss (function of `DFT`
        and `settlement_factor`, ref. [2], section 6.7.5.1.2, eqs. 27 and 28) and
        by the plastic loss (`PPL * Fv`, ref. [2], eq. 24).
        '''
        return _preload_after_losses(self)


    @cached_property
    def shell_force_at_rest (self):
        '''Shell force when no external loads are applied (Z1).

        The shell loads at rest are normally the self-weights of
        the structure supported by the flange.
        '''
        return self.Zg


    @cached_property
    def bolt_force_at_rest (self):
        '''Bolt axial force when no external loads are applied (Fs1).

        The bolt force at rest is just the bolt pretension.
        '''
        return self.bolt_preload


    @cached_property
    def bolt_moment_at_rest (self):
        '''Bolt bending moment when no external loads are applied (Ms1).'''

        Z1 = self.shell_force_at_rest
        Fs1 = self.bolt_force_at_rest
        return self._bolt_moment(Z1, Fs1)


    @cached_property
    def shell_force_at_small_displacement (self):
        '''Intermediate tower shell force (Z3).

        This auxiliary point gives the polynomial the right initial slope,
        evaluated according to ref. [2], Table G.1, point 3.
        '''

        # Retrieve point 2A (called point 0 in [1]) and determine Z as
        # a low pecentage of the theoretical pull Z0.
        Z0 = self._ideal_shell_force_at_tensile_ULS
        return 0.05 * self._stiffness_correction_factor * Z0


    @cached_property
    def bolt_force_at_small_displacement (self):
        '''Intermediate bolt axial force (Fs3).

        This auxiliary point gives the polynomial the right initial slope,
        evaluated according to ref. [2], Table G.1, point 3.
        '''

        # The slope between points P1 and P3 should match the
        # theoretical value of stiffness of the system.
        Z = self.shell_force_at_small_displacement
        return self.bolt_preload + self._polynomial_initial_slope * (Z - self.Zg)


    @cached_property
    def bolt_moment_at_small_displacement (self):
        '''Intermediate bolt bending moment (Ms3).'''

        Z3 = self.shell_force_at_small_displacement
        Fs3 = self.bolt_force_at_small_displacement
        return self._bolt_moment(Z3, Fs3)


    @cached_property
    def shell_force_at_tensile_ULS (self):
        '''Ultimate Limit State tower shell force (Z2).

        This is the tower shell force that brings the flange segment to its
        tensile ULS, evaluated according to ref. [2], Table G.1, point 2.
        '''

        Z0 = self._ideal_shell_force_at_tensile_ULS
        Z2 = self._stiffness_correction_factor * max(
            Z0 + self._total_gap_neutralization_shell_force,
            0.2 * Z0)
        return min(Z0, Z2)


    @cached_property
    def bolt_force_at_tensile_ULS (self):
        '''Bolt axial force at tensile failure (Fs2).

        This is the bolt's tensile force at the shell's ULS, assuming
        failure mode B (ref. [2], eq. G.36), scaled by the gap shape factor.
        '''

        # Bolt force at tensile ULS for sinusoidal gap shape
        Fs2 = self._ideal_bolt_force_at_tensile_ULS

        # Scale Fs2 based on the gap shape
        Fs1 = self.bolt_force_at_rest
        return Fs1 + (Fs2 - Fs1) * self.gap.shape_factor


    @cached_property
    def bolt_moment_at_tensile_ULS (self):
        '''Bolt bending moment at tensile failure (Ms2).'''

        Z2 = self.shell_force_at_tensile_ULS
        Fs2 = self.bolt_force_at_tensile_ULS
        return self._bolt_moment(Z2, Fs2)


    @cached_property
    def shell_force_at_closed_gap (self):
        '''Shell force necessary to completely close the imperfection gap (Z4).

        Evaluated according to ref. [2], eq. G.56.
        '''
        return self.shell_force_at_rest + self._total_gap_neutralization_shell_force

    @cached_property
    def _a_red (self):
        '''Equivalent flange inner length a', according to ref. [2], eq. G.12.'''
        return self.b / (self._prying_lever_ratio - 1)

    @cached_property
    def _a_star (self):
        '''Parameter a* of the bolt moment curve, according to ref. [2], eq. G.60.'''
        return max(0.4, min((self.t / (self._a_red + self.b))**2 , 1)) * self._a_red

    @cached_property
    def _total_gap_neutralization_shell_force (self):
        '''Incremental force to close the imperfection gap (delta-Z_gap,total).

        Evaluated according to ref. [2], eq. G.39, with the sign convention of
        this module: a negative value is a force closing the gap.
        '''
        MAXIMUM_GAP_NEUTRALIZATION_SHELL_FORCE = -1000 # 1 kN
        return min(self._parallel_gap_neutralization_shell_force - self._tilt_neutralization_shell_force,
                   MAXIMUM_GAP_NEUTRALIZATION_SHELL_FORCE)


    @cached_property
    def _parallel_gap_neutralization_shell_force (self):
        '''Force necessary to close a parallel imperfection gap (ref. [2], eq. G.40).'''
        s_avg = self.s * (1 + self.s_ratio) / 2
        c = self.central_angle * (self.R - s_avg/2)
        return -0.5 * self._gap_stiffness * self.gap.height * c


    @cached_property
    def _tilt_neutralization_shell_force (self):
        '''Force that neutralizes the flange tilt (DZ_gap,incl).

        This solves equation G.66 of ref. [2] and converts the resulting
        moment into a shell force with equation G.70.
        '''

        # The equation to be solved is A*[inv(B)*C]*M = D,
        # where A and B are 2x2 matrices, C and D are 2x1 matrices
        # and M is the scalar to be determined

        # General parameters
        s_mean = self.s * (1+self.s_ratio)/2
        w = self.a + self.b + self.s/2
        Af = w * self.t
        If = w * self.t**3 / 12
        R_sh = self.R - s_mean/2
        R_fl = self.R - w/2
        nu = self.E / (2*self.G) - 1
        k = self.E * s_mean**3 / (12*(1-nu**2))
        n = (3*(1-nu**2))**0.25 / (R_sh * s_mean)**0.5

        # Matrix A
        a11 = 1 / (2*k*n**3)
        a12 = 1 / (2*k*n**2)
        a21 = a12
        a22 = 1 / (k*n)
        A = np.array([[a11, a12],
                    [a21, a22]])

        # Matrix B
        b11 = a11 + R_fl * R_sh / (self.E * Af)
        b12 = a12
        b21 = a12
        b22 = a22 + R_fl * R_sh / (self.E * If)
        B = np.array([[b11, b12],
                    [b21, b22]])

        # Matrix C
        c1 = 0
        c2 = R_fl**2 / (self.E * If)
        C = np.array([c1, c2])

        # Matrix D
        d1 = 0.001    # ux: dummy value, sice irrelevant
        d2 = self.tilt_angle
        D = np.array([d1, d2])

        # Reduce matrix equation to H*M = D
        H = A @ (np.linalg.inv(B) @ C)

        # The last equation of the system H*M = D is h2*M = d2 and
        # we can therefore determine M from it as d2/h2.
        M = d2 / H[1]

        # Return the tower shell force force corresponding to M (ref. [2], eq. G.70)
        c = self.central_angle * R_sh
        b_mean = self.b + self.s/2 - s_mean/2
        return -c * M / (s_mean/2 + b_mean)

    
    @cached_property
    def _ideal_bolt_force_at_tensile_ULS (self):
        '''Bolt axial force at tensile failure for a sinusoidal gap shape (ref. [2], eq. G.36).'''
        return max(self.bolt.ultimate_tensile_capacity(), 1.25*self.bolt_force_at_rest)


    @cached_property
    def _ideal_shell_force_at_tensile_ULS (self):
        '''Tower shell force at the theoretical state of full prying (Z0, ref. [2], eq. G.37).'''
        return self._ideal_bolt_force_at_tensile_ULS / self._prying_lever_ratio


    @cached_property
    def _prying_lever_ratio (self):
        '''Lever ratio (a' + b) / a', according to ref. [2], section G.1.3.'''
        return _tobinaga_lever_ratio(self.a, self.b, self.t)


    @cached_property
    def _k_shell (self):
        '''Initial shell stiffness k_shell,ini, according to ref. [2], eqs. G.43 and G.44.'''
        s_avg = (self.s + self.s_ratio * self.s) / 2    # Average shell thickness
        L_gap = self.R * self.gap.angle                 # Gap lenght measured on the outer radius of the shell
        k_fac = max(1.8, 1.3 + (8.0e-4 - 1.6e-7 * (self.R*1000)) * (L_gap*1000))    # ref. [2], eq. G.44
        try:
            k_shell = float(self.k_shell)
        except Exception:
            k_shell = self.k_shell(self) if callable(self.k_shell) else self.E * s_avg / (k_fac * L_gap)

        log_data(self, L_gap=L_gap, k_fac=k_fac)
        return k_shell


    @cached_property
    def _k_flange (self):
        '''Flange bending stiffness k_fl, according to ref. [2], eq. G.46.'''
        L_gap = self.R * self.gap.angle                 # Gap lenght measured on the outer radius of the shell
        w = self.a + self.b + self.s/2      # flange segment length
        A = w * self.t                      # flange segment longitudinal cross-section area
        I = w * self.t**3 / 12              # flange segment longitudinal corss-section second moment of area
        EI = self.E * I
        GA = self.G * A
        L2 = L_gap**2
        k_flange = 384 * EI * GA / (L2 * (GA*L2 + 48*EI))   # ref. [2], eq. G.46
        log_data(self, A_cf=A, I_cf=I, k_fl=k_flange)
        return k_flange


    @cached_property
    def _gap_stiffness (self):
        '''Stiffness of the design gap.

        Returns the gap stiffness as a spring constant per unit of flange-segment
        arc length, calculated according to ref. [2], eqs. G.41 and G.42. The factor
        f_tot of eq. G.42 is further reduced for gap angles smaller than 30 deg.
        '''

        # Stiffness correction factor
        f_tot = min(1.0 + 1.5 * self.gap.angle/(pi/2), 2.5) * min(self.gap.angle/(pi/6), 1)**2

        log_data(self, f_tot=f_tot)

        # Total gap stiffness according to ref. [2], eq. G.41
        return f_tot * (self._k_shell + self._k_flange)


    @cached_property
    def _stiffness_correction_factor (self):
        '''Stiffness correction factor (alpha-k).

        This factor modifies the segment force to account for the effect of
        the gap spring, evaluated according to ref. [2], eqs. G.49 to G.51.
        '''

        from math import pi

        # Retrieve point P0
        Z0 = self._ideal_shell_force_at_tensile_ULS
        Fs0 = self._ideal_bolt_force_at_tensile_ULS

        # Evaluate the displacement u in the ultimate prying state.
        s_avg = self.s * (1 + self.s_ratio) / 2
        c = self.central_angle * (self.R - s_avg/2)
        I = c * self.t**3 / 12
        u = (Z0 * self.b**2 / (3 * self.E * I) + (Fs0 - self.bolt_preload) / (2 * self._bolt_axial_stiffness * self._a_red)) * (self._a_red + self.b)   # ref. [2], eq. G.51

        # Evaluate the segment stiffness
        s_avg = (self.s + self.s_ratio * self.s) / 2
        c = self.central_angle * (self.R - s_avg/2)
        k_seg = Z0 / (u * c)

        log_data(self, u=u, k_seg=k_seg)

        # Return the stiffness correction factor, acc. ref. [2], eq. G.49
        return min(1 + self._gap_stiffness / k_seg,
                   4 * pi/3 / self.gap.angle)


    @cached_property
    def _polynomial_initial_slope (self):
        '''Initial slope of the Fs(Z) polynomial.

        This slope is calculated according to ref. [2], eqs. G.53 and G.55,
        and scaled by the gap shape factor.
        '''

        # Load factor of the tension spring
        Ks = self._bolt_axial_stiffness
        Kp = self._flange_axial_stiffness
        p = Ks / (Ks + Kp)

        # Initial slope correction factor
        scf = -self._total_gap_neutralization_shell_force / (0.5 * self.bolt_preload)

        # Maximum allowable value of p
        p_max = (self.bolt_force_at_tensile_ULS - self.bolt_force_at_rest) / (self.shell_force_at_tensile_ULS - self.shell_force_at_rest)

        # Initial slope
        return min(scf*p, p_max) * self.gap.shape_factor


    # ----------------------------------------------------------------------
    #  DEPRECATED ATTRIBUTES AND METHODS
    #
    #  The following methods are deprecated and kept only for backwards
    #  compatibility
    # ----------------------------------------------------------------------

    @property
    def gap_height (self):
        '''Deprecated. Use `gap.height` instead.'''
        return self.gap.height

    @property
    def gap_angle (self):
        '''Deprecated. Use `gap.angle` instead.'''
        return self.gap.angle

    @property
    def gap_shape_factor (self):
        '''Deprecated. Use `gap.shape_factor` instead.'''
        return self.gap.shape_factor



class PolynomialLFlangeSegmentED2 (PolynomialLFlangeSegment):
    '''An L-Flange segment modeled as it will be modified in IEC 61400-6 ED2.

    This class has the same parameters and interface as `PolynomialLFlangeSegment`,
    but it replaces the `a*` parameter, the parallel-gap closing force and the
    gap stiffness with the formulations that will be introduced in
    IEC 61400-6 ED2 (ref. [5], section 9.7, changes #1 to #3).
    '''

    @cached_property
    def _a_star (self):
        '''Parameter a* of the bolt moment curve, as modified in IEC 61400-6 ED2.

        The lower bound is 0.3*a' instead of 0.4*a' (ref. [5], section 9.7, change #3).
        '''
        return max(0.3, min((self.t / (self._a_red + self.b))**2 , 1)) * self._a_red


    @cached_property
    def _parallel_gap_neutralization_shell_force (self):
        '''Force necessary to close a parallel imperfection gap, as modified in IEC 61400-6 ED2.

        The AMD1 force is multiplied by the torsional correction factor f_k,gap
        (ref. [5], sections 9.3.4 and 9.7, change #1).
        '''

        # Calaculate ED2 correction factor
        d = self.bolt.nominal_diameter
        alfa_ref = 1 - 0.630 * 2*d / (1.25*self.b + self.b)
        It_ref = alfa_ref * (1.25*self.b+self.b) * (2*d)**3 / 3
        alfa = 1 - 0.630 * self.t / (self.a + self.b)
        It = alfa * (self.a + self.b) * self.t**3 / 3
        f_kgap = max(1.0 , min( self.a/self.b / (It/It_ref)**0.5 , 1.5))

        # Return ED2 gap closing force
        return f_kgap * super()._parallel_gap_neutralization_shell_force


    @cached_property
    def _gap_stiffness (self):
        '''Stiffness of the design gap.

        Returns the gap stiffness as a spring constant per unit of flange-segment
        arc length, as it will be modified in IEC 61400-6 ED2 (ref. [5], section 9.7,
        change #2): the factor f_tot only applies to the shell stiffness. As for
        AMD1, f_tot is further reduced for gap angles smaller than 30 deg.
        '''

        # Stiffness correction factor
        f_tot = min(1.4 + 0.6 * self.gap.angle/(pi/2), 2.0) * min(self.gap.angle/(pi/6), 1)**2

        log_data(self, f_tot=f_tot)

        # Total gap stiffness according to ref. [5], section 9.7
        return f_tot * self._k_shell + self._k_flange



@dataclass
class PolynomialTFlangeSegment (PolynomialFlangeSegment):
    '''A T-Flange implementation of a `PolynomialFlangeSegment`.

    This class models a T-flange segment using the polynomial model from
    ref. [1] and [2]. It defines the reference points for the polynomial
    functions inherited from the parent class. The segment contains a pair of
    bolts, one on each side of the shell.

    The transfer functions are available through the `bolt_axial_force` and
    `bolt_bending_moment` methods, which accept either a scalar or a NumPy
    array of shell pull forces. Differently from the L-flange, the bolt
    bending moment in the tensile domain is a cubic polynomial (ref. [2],
    eq. G.61).

    Attributes:
        a (float): Distance from inner flange face to bolt hole center.
        b (float): Distance from bolt hole center to shell center-line.
        s (float): Shell thickness.
        t (float): Flange thickness.
        R (float): Shell outer curvature radius.
        central_angle (float): Angle subtended by the flange segment arc.
        Zg (float): Load on the shell at rest (e.g., dead weight), per bolt pair.
            Negative for compression.
        bolt (Bolt): The bolt in the flange segment.
        Fv (float): Design preload value before losses.
        Do (float): Bolt hole diameter.
        washer (Washer): The washer used. None if no washer.
        nut (Nut): The nut used.
        gap (Gap): The flange gap.
        tilt_angle (float): Flange tilt angle in radians. It has no default value,
            because ref. [2], section G.2.4, requires the inclination of T-flanges
            to be always considered: pass 0.0 explicitly only if the taper of the
            flange ensures that the first point of contact is directly below the
            tower wall (ref. [2], section 6.7.5.3).
        E (float): Young's modulus of the flange material.
        G (float): Shear modulus of the flange material.
        s_ratio (float): Ratio of bottom shell thickness to `s`.
        r (float): Radius of the rounding between shell and flange.
        k_shell (float, callable, None): Shell stiffness [N/m/m]. Can be a value, a
            function taking the flange segment as argument and returning the
            stiffness, or `None` to use the simplified formula of ref. [2], eq. G.43.
        PPL (float): Percentage of bolt pretension loss due to plastic strain,
            expressed as a fraction of `Fv` (e.g. 0.03 for 3%). Defaults to 0.0.
        DFT (float): Maximum allowable dry film thickness of coating applied on the contact
            area between washers and flange [m]. Defaults to 0.0 m.
        settlement_factor (float): Factor reducing the amount of bolt settlement, to be
            taken as 1.0 for torque tightening and 0.5 for tension tightening. Defaults to 1.0.
    '''

    a: float        # distance between inner face of the flange and center of the bolt hole
    b: float        # distance between center of the bolt hole and center-line of the shell
    s: float        # shell thickness
    t: float        # flange thickness
    R: float        # shell outer curvature radius
    central_angle: float     # angle subtended by the flange segment

    Zg: float       # load applied to the flange segment shell at rest (normally dead weight
                    # of tower + RNA, divided by the number of bolt pairs). Negative if compression.

    bolt: Bolt      # Bolt object representing the flange segment bolt
    Fv: float       # applied bolt preload before pretension losses

    Do: float       # Bolt hole diameter
    washer: Washer  # Bolt washer
    nut: Nut        # Bolt nut

    gap: Gap        # Flange gap object

    tilt_angle: float       # angle of flange tilt (no default: see class docstring)

    E: float = 210e9    # Young modulus of the flange
    G: float = 80.77e9  # Shear modulus of the flange
    s_ratio: float = 1.0    # Ratio of bottom shell thickness over s. Default s_botom = s.
    r: float = 0.01         # Rounding between flange and shell

    k_shell: float = None   # optional initial shell stiffness in [N/m/m].

    PPL: float = 0.0      # Percentage of pretension losses due to plastic strain

    DFT: float = 0.0        # Maximum allowable dry film thickness of coatings applied in the
                            # contact area between washers and flange [m]

    settlement_factor: float = 1.0   # Factor reducing the amount of bolt settlement, to be
                                     # taken as 1.0 for torque tightening and 0.5 for tension tightening.


    def failure_mode (self, fy_sh, fy_fl, gamma_0=1.1):
        '''Determines the failure mode of the flange.

        The ultimate shell force is evaluated for each of the failure modes
        A, B, D and E of ref. [2] (eqs. G.24 to G.29), and the governing mode
        is the one with the lowest ultimate shell force.

        Args:
            fy_sh (float): Yield stress of the shell material. Not used, since the
                plastic hinge of T-flanges is in the flange (ref. [2], Figure G.6).
            fy_fl (float): Yield stress of the flange material.
            gamma_0 (float, optional): Material factor. Defaults to 1.1.

        Returns:
            tuple: A tuple containing:

                - str: The failure mode ("A", "B", "D", or "E").
                - tuple: The shell force limits for each failure mode (Zu_A, Zu_B, Zu_D, Zu_E).
        '''

        from scipy.optimize import fsolve

        fd_fl = fy_fl / gamma_0

        Dw = self.washer.outer_diameter if self.washer else self.nut.bearing_diameter

        # Reduced distances between plastic hinges (ref. [2], eqs. G.28 and G.29)
        b_red = self.b - self.s/2 - 0.8*self.r
        b_E = self.b - (Dw + self.Do)/4 - self.s/2 - 0.8*self.r

        # Segment width at the middle surface of the shell (ref. [2], eq. G.34)
        # and segment width reduced by the bolt hole (ref. [2], eq. G.7)
        c = self.central_angle * (self.R - self.s/2)
        c_hole = c - self.Do

        # Failure mode A (ref. [2], eq. G.24)
        F_tRd = self.bolt.ultimate_tensile_capacity()   # Bolt ultimate tensile capacity
        Zu_A = 2*F_tRd                                  # Ultimate tower shell force for failure mode A

        # Flange cross-section ultimate capacity, ref. [2], eqs. G.8, G.10 and G.11
        Vu_fl = fd_fl * c * self.t / 3**0.5                # Pure shear ultimate capacity
        Mu_fl = fd_fl * c * self.t**2 / 4                  # Pure bending ultimate capacity
        MVu_fl = lambda V: Mu_fl * (1 - (V / Vu_fl)**2)**0.5 if V < Vu_fl else 0    # Bending ultimate capacity, concurrent with shear

        # Each half of the T-flange segment carries half of the shell force Z.
        # The iterative calculation of ref. [2], eq. G.11 is therefore done
        # for the half segment, given the function Zu_half(Mu_pl3) returning
        # the half shell force at flange shear ULS for a given value of Mu_pl3,
        # and the result is doubled (factor 2 of eqs. G.25 to G.27).
        def find_Zu (Zu_half):
            return 2 * fsolve(lambda Z: Zu_half(MVu_fl(Z)) - Z, Zu_half(Mu_fl))[0]

        # Failure mode B (ref. [2], eq. G.25)
        a_red = self._a_red
        Zu_B = find_Zu(lambda Mu_pl3: (F_tRd * a_red + Mu_pl3) / (a_red + b_red))

        # Failure mode D (ref. [2], eqs. G.7, G.9 and G.26)
        Mu_pl2 = fd_fl * c_hole * self.t**2 / 4
        DMu_pl2 = F_tRd/2 * (self.Do/2 + Dw/2)/2
        Zu_D = find_Zu(lambda Mu_pl3: (Mu_pl2 + DMu_pl2 + Mu_pl3) / b_red)

        # Failure mode E (ref. [2], eq. G.27)
        Zu_E = find_Zu(lambda Mu_pl3: (Mu_fl + Mu_pl3) / b_E)

        # Determine the governing failure mode
        Zu_min = min(Zu_A, Zu_B, Zu_D, Zu_E)

        if Zu_min == Zu_A:
            return "A",(Zu_A, Zu_B, Zu_D, Zu_E)
        elif Zu_min == Zu_B:
            return "B",(Zu_A, Zu_B, Zu_D, Zu_E)
        elif Zu_min == Zu_D:
            return "D",(Zu_A, Zu_B, Zu_D, Zu_E)
        elif Zu_min == Zu_E:
            return "E",(Zu_A, Zu_B, Zu_D, Zu_E)


    def validate (self, fy_sh, fy_fl, gamma_0=1.1):
        '''Verifies the applicability of the polynomial model.

        Checks if this T-Flange Segment fails according to failure mode A
        (ref. [2], eq. G.31).

        Args:
            fy_sh (float): Yield stress of the shell material.
            fy_fl (float): Yield stress of the flange material.
            gamma_0 (float, optional): Material factor. Defaults to 1.1.

        Raises:
            ValueError: If the failure mode is not "A".
        '''

        failure_mode, failure_shellforces = self.failure_mode(fy_sh, fy_fl, gamma_0=gamma_0)

        if failure_mode != "A":
            raise ValueError(f"Failure mode {failure_mode} detected, while only failure mode A is supported.")


    @cached_property
    def _bolt_axial_stiffness (self):
        '''The axial stiffness of the bolt.'''
        tw = self.washer.thickness if self.washer else 0
        return self.bolt.axial_stiffness(2*self.t + 2*tw)


    @cached_property
    def _bolt_bending_stiffness (self):
        '''The bending stiffness of the bolt.'''
        tw = self.washer.thickness if self.washer else 0
        return self.bolt.bending_stiffness(2*self.t + 2*tw)


    @cached_property
    def _flange_axial_stiffness (self):
        '''The axial stiffness of the flange in the thickness direction.'''
        # Stiffnes of flange w.r.t. compression in thickness direction,
        # when no gap is present. Calculated according to ref. [3] and [4].
        from math import pi, inf
        Dw = self.washer.outer_diameter if self.washer else self.nut.bearing_diameter
        Do = self.Do
        h = self.t * 2
        A = pi * ((Dw + h/10)**2 - Do**2) / 4
        Kf = self.E * A / h
        Kw = self.washer.axial_stiffness if self.washer else inf
        return 1 / (1/Kf + 2 * 1/Kw)


    def _bolt_moment (self, Z, Fs):
        '''Calculates the bolt bending moment, according to ref. [2], eqs. G.61 to G.63.

        Args:
            Z (float or numpy.polynomial.Polynomial): Shell pull force.
            Fs (float or numpy.polynomial.Polynomial): Bolt axial force.

        Returns:
            float or numpy.polynomial.Polynomial: The bolt bending moment.
        '''
        s_avg = self.s * (1 + self.s_ratio) / 2
        c = self.central_angle * (self.R - s_avg/2)
        A_tg= c *self.t
        I_tg = c * self.t**3 / 12

        ak = self._stiffness_correction_factor

        phi_T_low=(Fs-self.bolt_preload)*( self.b**2/(2*self.E*I_tg) + 1/(0.85*self.G*A_tg) )
        phi_T_high=Z*( self.b**2/(ak*4*self.E*I_tg) + 1/(0.85*self.G*A_tg) )

        Z0 = self._ideal_shell_force_at_tensile_ULS
        Z_close=self.shell_force_at_closed_gap

        bolt_rotation=phi_T_low+(phi_T_high-phi_T_low)/(Z0-Z_close)*(Z-Z_close)

        log_data(self, a_star=self._a_star)

        return bolt_rotation * 2*self._bolt_bending_stiffness


    @cached_property
    def _tensile_moment_polynomial (self):
        '''Polynomial Ms(Z) for the tensile domain.

        The bolt moment of ref. [2], eq. G.61 is the product of the bolt
        rotation interpolation, linear in Z, and of the bolt force polynomial,
        quadratic in Z. It is therefore a cubic polynomial, which is evaluated
        exactly here, instead of fitting a quadratic polynomial through the
        points P1, P2 and P3.

        Returns:
            numpy.polynomial.polynomial.Polynomial: The polynomial for the tensile domain.
        '''

        from numpy.polynomial.polynomial import Polynomial

        Fs = self._tensile_force_polynomial
        Z = Polynomial((0, 1), domain=Fs.domain, window=Fs.window, symbol="Z")
        return self._bolt_moment(Z, Fs)


    @cached_property
    def bolt_preload (self):
        '''float: Bolt preload after settlement and plastic losses.

        The preload `Fv` is reduced by the settlement loss (function of `DFT`
        and `settlement_factor`, ref. [2], section 6.7.5.1.2, eqs. 27 and 28) and
        by the plastic loss (`PPL * Fv`, ref. [2], eq. 24).
        '''
        return _preload_after_losses(self)


    @cached_property
    def shell_force_at_rest (self):
        '''Shell force when no external loads are applied (Z1).

        The shell loads at rest are normally the self-weights of
        the structure supported by the flange.
        '''
        return self.Zg

    @cached_property
    def bolt_force_at_rest (self):
        '''Bolt axial force when no external loads are applied (Fs1).

        The bolt force at rest is just the bolt pretension.
        '''
        return self.bolt_preload

    @cached_property
    def bolt_moment_at_rest (self):
        '''Bolt bending moment when no external loads are applied (Ms1).'''
        Z1 = self.shell_force_at_rest
        Fs1 = self.bolt_force_at_rest
        return self._bolt_moment(Z1, Fs1)

    @cached_property
    def shell_force_at_small_displacement (self):
        '''Intermediate tower shell force (Z3).

        This auxiliary point gives the polynomial the right initial slope,
        evaluated according to ref. [2], Table G.1, point 3.
        '''
        # Retrieve point 2A (called point 0 in [1]) and determine Z as
        # a low pecentage of the theoretical pull Z0.
        Z0 = self._ideal_shell_force_at_tensile_ULS
        return 0.05 * self._stiffness_correction_factor * Z0

    @cached_property
    def bolt_force_at_small_displacement (self):
        '''Intermediate bolt axial force (Fs3).

        This auxiliary point gives the polynomial the right initial slope,
        evaluated according to ref. [2], Table G.1, point 3.
        '''
        # The slope between points P1 and P3 should match the
        # theoretical value of stiffness of the system.
        Z = self.shell_force_at_small_displacement
        return self.bolt_preload + self._polynomial_initial_slope * (Z - self.Zg)

    @cached_property
    def bolt_moment_at_small_displacement (self):
        '''Intermediate bolt bending moment (Ms3).'''
        Z3 = self.shell_force_at_small_displacement
        Fs3 = self.bolt_force_at_small_displacement
        return self._bolt_moment(Z3, Fs3)

    @cached_property
    def shell_force_at_tensile_ULS (self):
        '''Ultimate Limit State tower shell force (Z2).

        This is the tower shell force that brings the flange segment to its
        tensile ULS, evaluated according to ref. [2], Table G.1, point 2.
        '''
        Z0 = self._ideal_shell_force_at_tensile_ULS
        Z2 = self._stiffness_correction_factor * max(
            Z0 + self._total_gap_neutralization_shell_force,
            0.2 * Z0)
        return min(Z0, Z2)

    @cached_property
    def bolt_force_at_tensile_ULS (self):
        '''Bolt axial force at tensile failure (Fs2).

        This is the bolt's tensile force at the shell's ULS, assuming
        failure mode A (ref. [2], eq. G.36), scaled by the gap shape factor.
        '''
        # Bolt force at tensile ULS for sinusoidal gap shape
        Fs2 = self._ideal_bolt_force_at_tensile_ULS

        # Scale Fs2 based on the gap shape
        Fs1 = self.bolt_force_at_rest
        return Fs1 + (Fs2 - Fs1) * self.gap.shape_factor

    @cached_property
    def bolt_moment_at_tensile_ULS (self):
        '''Bolt bending moment at tensile failure (Ms2).'''
        Z2 = self.shell_force_at_tensile_ULS
        Fs2 = self.bolt_force_at_tensile_ULS
        return self._bolt_moment(Z2, Fs2)

    @cached_property
    def shell_force_at_closed_gap (self):
        '''Shell force necessary to completely close the imperfection gap (Z4).

        Evaluated according to ref. [2], eq. G.56.
        '''
        return self.shell_force_at_rest + self._total_gap_neutralization_shell_force

    @cached_property
    def _a_red (self):
        '''Equivalent flange inner length a', according to ref. [2], eq. G.12.'''
        return self.b / (self._prying_lever_ratio - 1)

    @cached_property
    def _a_star (self):
        '''Parameter a* (ref. [2], eq. G.60). Not used by the T-flange model, only logged.'''
        return max(0.4, min((self.t / (self._a_red + self.b))**2 , 1)) * self._a_red

    @cached_property
    def _total_gap_neutralization_shell_force (self):
        '''Total incremental force to close the imperfection gap.

        Evaluated according to ref. [2], eq. G.39, with the sign convention of
        this module: a negative value is a force closing the gap.
        '''
        delta_Z_gap = self._parallel_gap_neutralization_shell_force
        delta_Z_gap_incl = self._tilt_neutralization_shell_force
        delta_Z_gap_c = self._early_prying_neutralization_shell_force

        return delta_Z_gap + delta_Z_gap_incl + delta_Z_gap_c

    @cached_property
    def _parallel_gap_neutralization_shell_force (self):
        '''Force necessary to close a parallel imperfection gap (ref. [2], eq. G.40).'''
        s_avg = self.s * (1 + self.s_ratio) / 2
        c = self.central_angle * (self.R - s_avg/2)
        return -0.5 * self._gap_stiffness * self.gap.height * c

    @cached_property
    def _tilt_neutralization_shell_force (self):
        '''Force that neutralizes the flange tilt (DZ_gap,incl).

        This solves equation G.66 of ref. [2] in closed form and converts the
        resulting moment into a shell force with equation G.69.
        '''
        # Calculate geometric properties
        Rm= self.R - self.s/2
        c_m=self.central_angle*Rm
        w=2*(self.a+self.b)

        # Shell parameters (ref. [2], eqs. G.67 and G.68)
        nu = self.E / (2*self.G) - 1
        K=self.E*self.s**3/(12*(1-nu**2))
        n=((3*(1-nu**2))**0.25)/((Rm*self.s)**0.5)

        # Second row of the matrix equation G.66, solved for M
        d11=1/(2*K*n**3)+Rm**2/(self.E*w*self.t)
        d22=1/(K*n)+12*Rm**2/(self.E*w*self.t**3)
        c12=1/(2*K*n**2)
        c22=1/(K*n)

        d02=self.tilt_angle * (d11*d22-c12**2)/(c22*d11-c12**2)
        M_tilt=d02*(self.E*w*self.t**3)/(12*Rm**2)
        delta_Z_gap_incl=-M_tilt*c_m/(w/2+self.b)

        log_data(self, M_tilt=M_tilt)

        return delta_Z_gap_incl

    @cached_property
    def _early_prying_neutralization_shell_force (self):
        '''Increase of the gap closing force due to early prying at the edges (DZ_gap,c).

        Evaluated according to ref. [2], eqs. G.47 and G.48.
        '''
        # Calculate geometric properties
        Rm= self.R - self.s/2
        c_m=self.central_angle*Rm
        A_tg= c_m *self.t
        I_tg = c_m * self.t**3 / 12

        # Initial linear gap stiffness, without the factor f_tot
        k_gap = self._k_shell + self._k_flange

        # Preload force at which edge contact occurs (ref. [2], eq. G.48)
        Fv_c=0.5*self.gap.height/( self.b**3/(3*self.E*I_tg) + self.b**2*self.a/(2*self.E*I_tg) + 1.2*self.b/(self.G*A_tg) + 1/(c_m*(k_gap)) )

        # Reduction due to early prying (ref. [2], eq. G.47)
        delta_Z_gap_c=-max((self.bolt_preload-Fv_c) * ((3*self.a+2*self.b)*self.b**2/(2*(self.a+self.b)**3)),0)

        log_data(self, Fv_c=Fv_c, A_tg=A_tg, I_tg=I_tg)

        return delta_Z_gap_c

    @cached_property
    def _ideal_bolt_force_at_tensile_ULS (self):
        '''Bolt axial force at tensile ULS (ref. [2], eq. G.36).'''
        # The bolt is loaded with the ultimate capacity, by definition.
        # If the pretension Fv is too close to the ultimate capacity Fsu, the polynomial
        # function may get too steep, therefore we make sure that Fs2 is at list 125% of Fv.
        return max(self.bolt.ultimate_tensile_capacity(), 1.25*self.bolt_force_at_rest)

    @cached_property
    def _ideal_shell_force_at_tensile_ULS (self):
        '''Shell force at the theoretical state of full separation (Z0, ref. [2], eq. G.38).'''
        return self._ideal_bolt_force_at_tensile_ULS * 2

    @cached_property
    def _prying_lever_ratio (self):
        '''Lever ratio (a' + b) / a', according to ref. [2], section G.1.3.'''
        return _tobinaga_lever_ratio(self.a, self.b, self.t)

    @cached_property
    def _k_shell (self):
        '''Initial shell stiffness k_shell,ini, according to ref. [2], eqs. G.43 and G.44.'''
        s_avg = (self.s + self.s_ratio * self.s) / 2    # Average shell thickness
        L_gap = self.R * self.gap.angle                 # Gap lenght measured on the outer radius of the shell
        k_fac = max(1.8, 1.3 + (8.0e-4 - 1.6e-7 * (self.R*1000)) * (L_gap*1000))    # ref. [2], eq. G.44
        try:
            k_shell = float(self.k_shell)
        except Exception:
            k_shell = self.k_shell(self) if callable(self.k_shell) else self.E * s_avg / (k_fac * L_gap)

        log_data(self, L_gap=L_gap, k_fac=k_fac, k_shell_ini=k_shell)
        return k_shell

    @cached_property
    def _k_flange (self):
        '''Flange bending stiffness k_fl, according to ref. [2], eq. G.46.'''
        L_gap = self.R * self.gap.angle     # Gap lenght measured on the outer radius of the shell
        w = (self.a + self.b)*2             # flange segment length
        A = w * self.t                      # flange segment longitudinal cross-section area
        I = w * self.t**3 / 12              # flange segment longitudinal corss-section second moment of area
        EI = self.E * I
        GA = self.G * A
        L2 = L_gap**2
        k_flange = 384 * EI * GA / (L2 * (GA*L2 + 48*EI))   # ref. [2], eq. G.46
        log_data(self, A_cf=A, I_cf=I, k_fl=k_flange)
        return k_flange

    @cached_property
    def _gap_stiffness (self):
        '''Stiffness of the design gap.

        Returns the gap stiffness as a spring constant per unit of flange-segment
        arc length, calculated according to ref. [2], eqs. G.41 and G.42.
        '''
        f_tot = min(1.0 + 1.5 * self.gap.angle/(pi/2), 2.5)

        log_data(self, f_tot=f_tot)

        # Total gap stiffness according to ref. [2], eq. G.41
        return f_tot * (self._k_shell + self._k_flange)

    @cached_property
    def _stiffness_correction_factor (self):
        '''Stiffness correction factor (alpha-k).

        This factor modifies the segment force to account for the effect of
        the gap spring, evaluated according to ref. [2], eqs. G.49, G.50 and G.52.
        '''

        from math import pi

        Rm= self.R - self.s/2
        cm=self.central_angle*Rm

        # Retrieve point 0
        Z0 = self._ideal_shell_force_at_tensile_ULS
        Fs0 = self._ideal_bolt_force_at_tensile_ULS

        # Evaluate the displacement u in the ultimate prying state (ref. [2], eq. G.52)
        u = (Fs0 - self.bolt_preload) / (2*self._bolt_axial_stiffness)

        # Evaluate the segment stiffness (ref. [2], eq. G.50)
        k_seg = Z0 / (u * cm)

        log_data(self, u=u, k_seg=k_seg)

        # Return the stiffness correction factor, acc. ref. [2], eq. G.49
        return min(1 + self._gap_stiffness / k_seg,
                   4 * pi/3 / self.gap.angle)

    @cached_property
    def _polynomial_initial_slope (self):
        '''Initial slope of the Fs(Z) polynomial.

        This slope is calculated according to ref. [2], eqs. G.54 and G.55,
        and scaled by the gap shape factor.
        '''

        # Load factor of the tension spring
        p = self._bolt_axial_stiffness / (self._bolt_axial_stiffness + self._flange_axial_stiffness)

        # Maximum allowable value of the slope
        p_max = (self.bolt_force_at_tensile_ULS - self.bolt_force_at_rest) / (self.shell_force_at_tensile_ULS - self.shell_force_at_rest)

        # Gap closing ratio, considering the preload of both bolts of the segment
        gap_closing_ratio = -self._total_gap_neutralization_shell_force / (2 * self.bolt_preload * 0.5)

        # Initial slope
        return min(0.5 * p * gap_closing_ratio, p_max) * self.gap.shape_factor


    # ----------------------------------------------------------------------
    #  DEPRECATED ATTRIBUTES AND METHODS
    #
    #  The following methods are deprecated and kept only for backwards
    #  compatibility
    # ----------------------------------------------------------------------

    @property
    def gap_height (self):
        '''Deprecated. Use `gap.height` instead.'''
        return self.gap.height

    @property
    def gap_angle (self):
        '''Deprecated. Use `gap.angle` instead.'''
        return self.gap.angle

    @property
    def gap_shape_factor (self):
        '''Deprecated. Use `gap.shape_factor` instead.'''
        return self.gap.shape_factor



class PolynomialTFlangeSegmentED2 (PolynomialTFlangeSegment):
    '''A T-Flange segment modeled as it will be modified in IEC 61400-6 ED2.

    This class has the same parameters and interface as `PolynomialTFlangeSegment`,
    but it replaces the gap stiffness with the formulation that will be introduced
    in IEC 61400-6 ED2 (ref. [5], sections 9.3.3 and 9.7, change #2). The other
    ED2 changes of ref. [5], section 9.7 (torsional correction factor and a*)
    only concern L-flanges.
    '''

    @cached_property
    def _gap_stiffness (self):
        '''Stiffness of the design gap.

        Returns the gap stiffness as a spring constant per unit of flange-segment
        arc length, as it will be modified in IEC 61400-6 ED2 (ref. [5], section 9.7,
        change #2): the factor f_tot only applies to the shell stiffness.
        '''
        f_tot = min(1.4 + 0.6 * self.gap.angle/(pi/2), 2.0)

        log_data(self, f_tot=f_tot)

        # Total gap stiffness according to ref. [5], section 9.7
        return f_tot * self._k_shell + self._k_flange




def _tobinaga_lever_ratio (a, b, t):
    '''Lever ratio (a' + b) / a' of a flange segment, according to ref. [2], section G.1.3.

    The equivalent flange inner length a' is evaluated with the action point
    correction of Tobinaga & Ishihara (ref. [6]; ref. [2], eqs. G.12 to G.15).
    The correction is only applied within the validity range of ref. [2],
    eqs. G.16 and G.17: outside that range a' = a is used, and a warning is
    logged when 1.25 < a/b.

    Args:
        a (float): Distance from inner flange face to bolt hole center.
        b (float): Distance from bolt hole center to shell center-line.
        t (float): Flange thickness.

    Returns:
        float: The lever ratio (a' + b) / a'.
    '''

    rigid_lever_ratio = (a + b) / a

    # The flange is considered rigid if a/b <= 1.25: the reaction force acts
    # at the flange edge, therefore a' = a.
    if a/b <= 1.25:
        return rigid_lever_ratio

    # Aspect ratio of the flange section (ref. [2], eq. G.14)
    alpha = t / (a + b)

    # Out of the validity range of eqs. G.16 and G.17, a' = a
    if a/b > 2.25:
        logger.warning("Flange geometry out of validity range for Tobinaga & Ishihara action point correction (a/b > 2.25): a' = a is used.")
        return rigid_lever_ratio
    if not (0.55 - 0.12 * (a/b) <= alpha <= 1.0):
        logger.warning("Flange geometry out of validity range for Tobinaga & Ishihara action point correction (-0.12*a/b + 0.55 <= t/(a+b) <= 1 not satisfied): a' = a is used.")
        return rigid_lever_ratio

    # Modification factor of a (ref. [2], eqs. G.13 and G.15)
    beta = (a/b - 1.25)**0.32 + 0.45
    lambda_ = 1 - (1 - alpha**beta)**5

    # Equivalent flange inner length (ref. [2], eq. G.12)
    a_red = lambda_ * a

    return (a_red + b) / a_red



def _preload_after_losses (fseg):
    '''Bolt preload of a flange segment after settlement and plastic losses.

    The preload `fseg.Fv` is reduced by the settlement loss (ref. [2],
    section 6.7.5.1.2, eqs. 27 and 28) and by the plastic loss `fseg.PPL * fseg.Fv`
    (ref. [2], eq. 24).

    Args:
        fseg (FlangeSegment): A flange segment with the attributes `Fv`, `PPL`,
            `DFT` and `settlement_factor`, and the properties `_bolt_axial_stiffness`
            and `_flange_axial_stiffness`.

    Returns:
        float: The preload after losses.
    '''
    from metrum.units import um, kN

    # Preload losses due to settlements
    if fseg.DFT > 0:
        DFT_ref = 1*um
        f_z_tot = max( (8*(fseg.DFT/DFT_ref)**(-0.6))*fseg.DFT , 50*um )
        bolt_resilience = 1 / fseg._bolt_axial_stiffness
        flange_resilience = 1 / fseg._flange_axial_stiffness
        delta_Fv_z = fseg.settlement_factor * (f_z_tot / (bolt_resilience + flange_resilience))
    else:
        delta_Fv_z = 0*kN

    # Preload losses due to plastic deformation
    delta_Fv_pl = fseg.PPL * fseg.Fv

    return fseg.Fv - delta_Fv_z - delta_Fv_pl



def bolt_markov_matrix (fseg, flange_markov_matrix, bending_factor=0.0, macro_geometric_factor=1.0, mean_factor=1.0, range_factor=1.0, SMF=1.0):
    '''Evaluates the bolt Markov matrix from the shell Markov matrix.

    The flange bending moment ranges and means are converted into shell pull
    forces acting on the flange segment, and then into minimum and maximum bolt
    stresses through the transfer functions of the flange segment. The bolt
    stress is the sum of the axial stress and of the bending stress scaled by
    `bending_factor`, evaluated on the bolt thread cross-section, according
    to ref. [2], section 6.7.5.4, eq. (42).

    The flange bending moment M is converted into the segment shell force
    as Z = M / W * A_seg (ref. [2], eq. G.33), where W is the section modulus
    of the tower shell, evaluated at the shell mid-surface, and A_seg the
    shell cross-section area of the segment.

    Args:
        fseg (PolynomialFlangeSegment): The flange segment containing the
            transfer functions. **WARNING:** Currently only works with
            PolynomialFlangeSegments.
        flange_markov_matrix (pyflange.fatigue.MarkovMatrix): The Markov
            matrix of bending moments on the flange.
        bending_factor (float, optional): Factor for the bending portion of the total
            stress range (kappa_b of ref. [2], eq. (43)). Defaults to 0.0
            (axial stress only).
        macro_geometric_factor (float, optional): Factor for macro geometric influences
            (the SCF of ref. [2], section 6.7.5.3), affecting dead weight, mean,
            and range values. If omitted, it defaults to 1.0.
        mean_factor (float, optional): Factor for the mean values of the tower's
            bending moments. If omitted, it defaults to 1.0.
        range_factor (float, optional): Factor for the range of the tower's bending
            moments. If omitted, it defaults to 1.0.
        SMF (float, optional): Stress Multiplication Factor to be applied 
            to bolt Markov matrix stress ranges. If omitted, it defaults to 1.0.

    Returns:
        pyflange.fatigue.MarkovMatrix: The bolt's Markov matrix representing
            the stress history. Its cycles and duration are the same as the
            flange Markov matrix.
    '''
    from .fatigue import MarkovMatrix

    #Flange Geometry
    Rm = fseg.R - fseg.s/2
    flange_W = pi/4 * (fseg.R**4 - (fseg.R-fseg.s)**4) / Rm
    shell_A = fseg.s * fseg.central_angle * Rm

    # Bolt Geometry
    bolt_A = fseg.bolt.thread_cross_section.area
    bolt_W = fseg.bolt.thread_cross_section.elastic_section_modulus

    # Shell Markov Matrix
    Z_cycles = flange_markov_matrix.cycles
    Z_mean   = flange_markov_matrix.mean / flange_W * shell_A * mean_factor * macro_geometric_factor
    Z_ranges = flange_markov_matrix.range / flange_W * shell_A * range_factor * macro_geometric_factor
    Z_min = Z_mean - Z_ranges/2 + fseg.Zg * macro_geometric_factor
    Z_max = Z_mean + Z_ranges/2 + fseg.Zg * macro_geometric_factor

    # Bolt ess Markov Matrix
    S_min = fseg.bolt_axial_force(Z_min)/bolt_A + bending_factor * fseg.bolt_bending_moment(Z_min)/bolt_W
    S_max = fseg.bolt_axial_force(Z_max)/bolt_A + bending_factor * fseg.bolt_bending_moment(Z_max)/bolt_W
    return MarkovMatrix(range = SMF * np.abs(S_max - S_min),
                        mean = (S_min + S_max) / 2,
                        cycles = Z_cycles,
                        duration = flange_markov_matrix.duration)



def shell_stiffness (shell_radius, shell_thickness, gap_angle):
    '''Calculates the axial stiffness of a flange-segment shell via interpolation.

    This function interpolates shell stiffness from a dataset of 2160 FEA results,
    covering 8 diameters (3 m to 10 m), 15 shell thicknesses (10 mm to 150 mm)
    and 18 gap angles (10° to 180°). This is the database of pre-computed
    stiffness values described in ref. [5], section 9.3.1, as an alternative
    to the simplified formula of ref. [2], eqs. G.43 and G.44.
    Values outside the range covered by the dataset are not extrapolated, and
    `nan` is returned instead.

    Example:
        ```python
        from math import pi
        from pyflange.flangesegments import shell_stiffness

        k_shell = shell_stiffness(4.0, 0.054, pi/6)   # R = 4 m, s = 54 mm, 30° gap
        ```

    Args:
        shell_radius (float): Outer radius of the shell.
        shell_thickness (float): Thickness of the shell.
        gap_angle (float): Angular length of the flange gap in radians.

    Returns:
        float: The interpolated shell stiffness (K_sh), or `nan` if the
            arguments are outside the dataset range.
    '''
    interpolate_shell_stiffness = _load_shell_stiffness_interpolator()
    return interpolate_shell_stiffness(shell_radius, shell_thickness, gap_angle)


@cache
def _load_shell_stiffness_interpolator ():
    # Loads and caches the shell stiffness interpolator.
    # Returns scipy.interpolate.LinearNDInterpolator: The interpolator function.

    from .utils import load_csv_database
    from scipy.interpolate import LinearNDInterpolator
    ssdf = load_csv_database("flangesegments.shell_stiffness", index_col=None)
    return LinearNDInterpolator(list(zip(ssdf['shell_radius'], ssdf['shell_thickness'], ssdf['gap_angle']/180*pi)), ssdf['shell_stiffness'])
