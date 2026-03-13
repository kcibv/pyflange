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

"""Fatigue calculation tools.

This module defines functions and classes to support structural fatigue
calculations for large flanges. It includes utilities for reading Markov matrices,
defining various types of S-N curves (fatigue curves), and performing
bolt fatigue analysis.
"""

import numpy as np
import pandas as pd
from math import sqrt, pi, tan, exp, log10

from dataclasses import dataclass
import functools

from .flangesegments import FlangeSegment


@dataclass
class MarkovMatrix:
    """A Markov Matrix representing a load history.

    This class represents a load history in the form of a Markov Matrix.
    The 'load' can be expressed in terms of moments, forces, stresses, etc.

    Attributes:
        range (np.ndarray): Array of stress ranges (or load ranges in general).
        mean (np.ndarray): Array of mean stress values (or load values in general).
        cycles (np.ndarray): Array of cycles corresponding to each load value.
        duration (float): Duration of the load history. Defaults to 1.
    """

    range: np.ndarray
    mean: np.ndarray
    cycles: np.ndarray
    duration: float = 1


class FatigueCurve:
    """A Wöhler (S-N) curve base class.

    This is a base class for creating Wohler curves. It should not be
    instantiated directly.
    """

    def N(self, DS):
        """Calculates the number of cycles to failure for a given stress range.

        Args:
            DS (float or np.ndarray): The stress range(s).

        Returns:
            float or np.ndarray: The number of cycles that produce a fatigue failure.
        """
        pass

    def DS(self, N):
        """Calculates the stress range for a given number of cycles.

        Args:
            N (float or np.ndarray): The number of cycles to failure.

        Returns:
            float or np.ndarray: The corresponding stress range that produces
                fatigue failure.
        """
        pass

    def damage(self, n, DS):
        """Calculates the fatigue damage for a given number of cycles and stress range.

        Calculates the damage as D = n / N(DS).

        Args:
            n (float or np.ndarray): The number of applied cycles.
            DS (float or np.ndarray): The stress range(s).

        Returns:
            float or np.ndarray: The fatigue damage.
        """
        return n / self.N(DS)

    def cumulated_damage(self, markov_matrix):
        """Calculates the cumulated damage according to Miner's rule.

        Args:
            markov_matrix (MarkovMatrix): The load history expressed as a
                MarkovMatrix object.

        Returns:
            float: The cumulated fatigue damage under the given load history.
        """

        n = markov_matrix.cycles    # array of number of cycles
        DS = markov_matrix.range    # array of stress ranges
        D = self.damage(n, DS)      # array of damages
        return np.nansum(D)         # total damage


@dataclass
class SingleSlopeFatigueCurve(FatigueCurve):
    """Wöhler curve with a single logarithmic slope.

    This class implements the FatigueCurve interface for a curve with a single
    slope m.

    Attributes:
        m (float): The logarithmic slope of the fatigue curve.
        DS_ref (float): Arbitrary reference stress range.
        N_ref (float): The number of cycles that produce failure under the
            stress range DS_ref.
    """

    m: float
    DS_ref: float
    N_ref: float

    @functools.cached_property
    def a(self):
        """float: The constant 'a' in the S-N curve equation N * DS^m = a."""
        return self.DS_ref ** self.m * self.N_ref

    def N(self, DS):
        return self.a / DS**self.m

    def DS(self, N):
        return (self.a / N)**(1/self.m)


class MultiSlopeFatigueCurve(FatigueCurve):
    """Multi-Slope Fatigue Curve.

    This class represents a FatigueCurve composed of multiple slopes.
    It takes any number of SingleSlopeFatigueCurve objects and uses the
    maximum number of cycles (or stress range) among them for a given input.

    Attributes:
        curves (tuple): A tuple of SingleSlopeFatigueCurve objects.
    """

    def __init__(self, *fatigue_curves):
        """Initializes the MultiSlopeFatigueCurve.

        Args:
            *fatigue_curves (SingleSlopeFatigueCurve): Any number of fatigue curve
                components.
        """
        self.curves = fatigue_curves

    def N(self, DS):
        """Calculates the number of cycles to failure.

        Args:
            DS (float or np.ndarray): The stress range(s).

        Returns:
            float or np.ndarray: The maximum number of cycles from the component curves.
        """
        return np.maximum.reduce([curve.N(DS) for curve in self.curves])

    def DS(self, N):
        """Calculates the stress range.

        Args:
            N (float or np.ndarray): The number of cycles.

        Returns:
            float or np.ndarray: The maximum stress range from the component curves.
        """
        return np.maximum.reduce([curve.DS(N) for curve in self.curves])


class DoubleSlopeFatigueCurve(MultiSlopeFatigueCurve):
    """Wöhler curve with double logarithmic slope.

    This class implements the FatigueCurve interface for a curve with two
    slopes m1 and m2, intersecting at a specific point (DS12, N12).
    """

    def __init__(self, m1, m2, DS12, N12):
        """Initializes the DoubleSlopeFatigueCurve.

        Args:
            m1 (float): The logarithmic slope of the lower cycle values.
            m2 (float): The logarithmic slope of the higher cycle values.
            DS12 (float): The stress range at the knee point where slopes meet.
            N12 (float): The number of cycles at the knee point where slopes meet.
        """
        curve1 = SingleSlopeFatigueCurve(m1, DS12, N12)
        curve2 = SingleSlopeFatigueCurve(m2, DS12, N12)
        super().__init__(curve1, curve2)


class BoltFatigueCurve(DoubleSlopeFatigueCurve):
    """Bolt Fatigue Curve according to IEC 61400-6 AMD1.

    Creates a DoubleSlopeFatigueCurve with logarithmic slopes m1=3 and m2=5,
    a knee point at 2 million cycles, and a stress range depending on the
    bolt diameter as specified by IEC 61400-6 AMD1.

    Args:
        diameter (float): The bolt diameter in meters.
        DS_ref (float, optional): Reference stress range at 2e6 cycles.
            Defaults to 50 MPa (50e6 Pa).
        gamma_M (float, optional): The material factor. Defaults to 1.1.
    """

    def __init__(self, diameter, DS_ref=50e6, gamma_M=1.1):
        N12 = 2.0e6    # knee point
        m1 = 3
        m2 = 5
        if diameter <= 0.030:
            DSc = DS_ref    # reference stress range, in Pa
        elif diameter <= 0.072:
            DSc = DS_ref * (0.030/diameter)**0.1
        else:
            DSc = DS_ref * (0.030/diameter)**0.1 * (0.072/diameter)**0.25

        # Delegate the rest of the initialization to the parent class
        super().__init__(m1, m2, DSc/gamma_M, N12)


@dataclass
class BoltFatigueAnalysis:
    """Fatigue analysis for a flange bolt.

    Calculates fatigue damage and life for a bolt given its flange segment model,
    load history, and fatigue curve.

    Attributes:
        fseg (FlangeSegment): The FlangeSegment object for the analysis.
        flange_mkvm (MarkovMatrix): The Markov matrix representing the bending
            moments history on the flange.
        custom_fatigue_curve (FatigueCurve, optional): The fatigue curve for the
            bolt. If None, a BoltFatigueCurve is generated based on the bolt
            diameter. Defaults to None.
        allowable_damage (float): The maximum allowable fatigue damage.
            Defaults to 1.0.
        SMF (float): Stress Multiplication Factor applied to bolt stress ranges.
            Defaults to 1.0.
    """

    fseg: FlangeSegment
    flange_mkvm: MarkovMatrix
    custom_fatigue_curve: FatigueCurve = None
    allowable_damage: float = 1.0
    SMF: float = 1.0

    @functools.cached_property
    def fatigue_curve(self):
        """FatigueCurve: The fatigue curve used in the analysis."""
        return self.custom_fatigue_curve or BoltFatigueCurve(self.fseg.bolt.nominal_diameter)

    @functools.cached_property
    def bolt_mkvm(self):
        """MarkovMatrix: The Markov matrix of axial+bending stress ranges in the bolt."""
        from math import log
        from .flangesegments import bolt_markov_matrix
        bending_factor = max(0.5, 0.5 + 0.5*log(self.fseg.bolt.nominal_diameter/0.036) / log(150/36))
        return bolt_markov_matrix(self.fseg, self.flange_mkvm, bending_factor, SMF=self.SMF)

    @functools.cached_property
    def damage(self):
        """float: The total cumulated fatigue damage."""
        return self.fatigue_curve.cumulated_damage(self.bolt_mkvm)

    @functools.cached_property
    def fatigue_life(self):
        """float: The fatigue life in the same units as flange_mkvm.duration."""
        return self.allowable_damage / self.damage * self.flange_mkvm.duration
