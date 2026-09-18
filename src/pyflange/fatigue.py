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

r"""Fatigue calculation and S-N curve tools for flanged connections.

This module provides data structures, fatigue (Wöhler) curve representations,
and analysis classes to evaluate structural fatigue damage and fatigue life
for large bolted flanges.

Key components provided by this module:

- `MarkovMatrix`: Encapsulates a discrete cyclic load or stress history (load
  ranges, mean values, cycle counts, and duration) obtained from rainflow counting
  or design load spectra.
- `FatigueCurve`: Abstract base class defining the S-N curve interface and
  Palmgren-Miner linear cumulative damage accumulation:
  $$D = \sum_{i} \frac{n_i}{N(\Delta S_i)}$$
- `SingleSlopeFatigueCurve`: Represents a Wöhler curve with a single logarithmic
  slope exponent $m$, defined by the power law $N \cdot \Delta S^m = a$.
- `MultiSlopeFatigueCurve`: Represents a composite multi-slope S-N curve
  constructed as the upper envelope of multiple single-slope curve segments.
- `DoubleSlopeFatigueCurve`: Represents a bilinear S-N curve with two slopes
  ($m_1, m_2$) intersecting at a transition knee point $(\Delta S_{12}, N_{12})$.
- `BoltFatigueCurve`: Specialized bilinear S-N curve for threaded bolts conforming
  to IEC 61400-6 AMD1, incorporating nominal diameter size-effect corrections
  and a material partial safety factor $\gamma_M$.
- `BoltFatigueAnalysis`: End-to-end fatigue analysis pipeline that evaluates bolt
  stress spectra from flange bending moment histories, calculates total cumulative
  damage, and predicts bolt fatigue life.

Units:
    All parameters must be supplied in a consistent system of units. When using
    SI units:
    - Distances and diameters: meters (m)
    - Forces: newtons (N)
    - Bending moments: newton-meters (N·m)
    - Stresses and pressures: pascals (Pa)
    - Duration and fatigue life: consistent time units (e.g., seconds or years)
"""

from __future__ import annotations

from dataclasses import dataclass
import functools
import numpy as np

from .flangesegments import FlangeSegment


@dataclass
class MarkovMatrix:
    r"""A Markov matrix representation of a cyclic load history.

    Encapsulates a discrete load or stress spectrum resulting from cycle counting
    (e.g., rainflow counting). The load history is characterized by arrays of cycle
    ranges, mean values, and associated cycle counts, along with an overall duration.

    The physical quantity represented can be stress (Pa), force (N), bending
    moment (N·m), etc., provided that consistent units are maintained throughout
    downstream calculations.

    Attributes:
        range (np.ndarray): Array of cycle load or stress ranges ($\Delta S$ or $\Delta L$).
        mean (np.ndarray): Array of mean load or stress values ($S_m$ or $L_m$).
        cycles (np.ndarray): Array containing the number of applied cycles ($n$) for
            each (range, mean) bin.
        duration (float): Time duration represented by the load history (e.g., in
            seconds or years). Defaults to 1.0.
    """

    range: np.ndarray
    mean: np.ndarray
    cycles: np.ndarray
    duration: float = 1.0


class FatigueCurve:
    r"""Abstract base class for Wöhler (S-N) fatigue strength curves.

    This class provides the interface for evaluating fatigue endurance (number of
    cycles to failure $N$ as a function of stress range $\Delta S$, and vice versa),
    as well as calculating fatigue damage and cumulative damage using the
    Palmgren-Miner linear cumulative damage hypothesis:

    $$D = \sum_{i} \frac{n_i}{N(\Delta S_i)}$$

    This class should not be instantiated directly; use one of its concrete
    subclasses such as `SingleSlopeFatigueCurve`, `DoubleSlopeFatigueCurve`, or
    `BoltFatigueCurve`.
    """

    def N(self, DS: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the number of cycles to failure for a given stress range.

        Args:
            DS (float | np.ndarray): Applied stress range(s) $\Delta S$.

        Returns:
            float | np.ndarray: Number of cycles to fatigue failure $N(\Delta S)$.
        """
        pass

    def DS(self, N: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the allowable stress range for a given number of cycles.

        Args:
            N (float | np.ndarray): Number of cycles to failure $N$.

        Returns:
            float | np.ndarray: Allowable stress range(s) $\Delta S(N)$ producing
                fatigue failure at $N$ cycles.
        """
        pass

    def damage(self, n: float | np.ndarray, DS: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates fatigue damage according to the Palmgren-Miner rule.

        Computes the elementary fatigue damage ratio:

        $$D = \frac{n}{N(\Delta S)}$$

        Args:
            n (float | np.ndarray): Number of applied load cycles.
            DS (float | np.ndarray): Applied stress range(s) $\Delta S$.

        Returns:
            float | np.ndarray: The resulting fatigue damage ratio(s).
        """
        return n / self.N(DS)

    def cumulated_damage(self, markov_matrix: MarkovMatrix) -> float:
        r"""Calculates cumulative fatigue damage according to Miner's sum.

        Evaluates the total cumulative damage over all bins in the load spectrum:

        $$D_{\text{total}} = \sum_{i} \frac{n_i}{N(\Delta S_i)}$$

        Bins producing `NaN` damage values (e.g., zero stress or invalid evaluations)
        are ignored in the summation via `numpy.nansum`.

        Args:
            markov_matrix (MarkovMatrix): Load history represented as a Markov matrix
                containing arrays of cycle counts and stress ranges.

        Returns:
            float: Total cumulative fatigue damage index $D_{\text{total}}$.
        """
        n = markov_matrix.cycles    # array of number of cycles
        DS = markov_matrix.range    # array of stress ranges
        D = self.damage(n, DS)      # array of damages
        return float(np.nansum(D))  # total damage


@dataclass
class SingleSlopeFatigueCurve(FatigueCurve):
    r"""Wöhler (S-N) fatigue curve with a single logarithmic slope.

    Represents a linear S-N relationship in log-log space governed by the power-law
    equation:

    $$N \cdot \Delta S^m = a = \Delta S_{\text{ref}}^m \cdot N_{\text{ref}}$$

    where $m$ is the inverse slope exponent, $\Delta S_{\text{ref}}$ is a reference
    stress range, and $N_{\text{ref}}$ is the corresponding cycle count at failure.

    Attributes:
        m (float): Negative inverse slope of the S-N curve in log-log coordinates
            ($-\frac{\mathrm{d}\log N}{\mathrm{d}\log \Delta S}$).
        DS_ref (float): Reference stress range $\Delta S_{\text{ref}}$ (e.g., in Pa).
        N_ref (float): Reference number of cycles to failure $N_{\text{ref}}$ at
            stress range `DS_ref`.
    """

    m: float
    DS_ref: float
    N_ref: float

    @functools.cached_property
    def a(self) -> float:
        r"""float: Curve capacity constant $a = \Delta S_{\text{ref}}^m \cdot N_{\text{ref}}$."""
        return self.DS_ref ** self.m * self.N_ref

    def N(self, DS: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the number of cycles to failure for given stress range(s).

        Evaluates:

        $$N = \frac{a}{\Delta S^m}$$

        Args:
            DS (float | np.ndarray): Applied stress range(s) $\Delta S$.

        Returns:
            float | np.ndarray: Number of cycles to failure $N$.
        """
        return self.a / DS**self.m

    def DS(self, N: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the allowable stress range for given cycle count(s).

        Evaluates:

        $$\Delta S = \left(\frac{a}{N}\right)^{1/m}$$

        Args:
            N (float | np.ndarray): Number of cycles to failure $N$.

        Returns:
            float | np.ndarray: Allowable stress range(s) $\Delta S$.
        """
        return (self.a / N)**(1/self.m)


class MultiSlopeFatigueCurve(FatigueCurve):
    r"""Multi-slope composite Wöhler (S-N) fatigue curve.

    Represents an S-N curve formed by combining multiple `SingleSlopeFatigueCurve`
    segments. For any given stress range or cycle count, the effective capacity is
    evaluated as the upper envelope (maximum endurance / maximum allowable stress
    range) among all constituent curve segments.

    Attributes:
        curves (tuple[SingleSlopeFatigueCurve, ...]): Tuple of constituent single-slope
            fatigue curve segments.
    """

    def __init__(self, *fatigue_curves: SingleSlopeFatigueCurve) -> None:
        """Initializes the MultiSlopeFatigueCurve.

        Args:
            *fatigue_curves (SingleSlopeFatigueCurve): One or more single-slope
                fatigue curve segments composing the multi-slope curve.
        """
        self.curves = fatigue_curves

    def N(self, DS: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the maximum number of cycles to failure among all curve segments.

        Evaluates:

        $$N(\Delta S) = \max_{i} N_i(\Delta S)$$

        Args:
            DS (float | np.ndarray): Applied stress range(s) $\Delta S$.

        Returns:
            float | np.ndarray: Maximum number of cycles to failure across all component
                curves.
        """
        return np.maximum.reduce([curve.N(DS) for curve in self.curves])

    def DS(self, N: float | np.ndarray) -> float | np.ndarray:
        r"""Calculates the maximum allowable stress range among all curve segments.

        Evaluates:

        $$\Delta S(N) = \max_{i} \Delta S_i(N)$$

        Args:
            N (float | np.ndarray): Number of cycles $N$.

        Returns:
            float | np.ndarray: Maximum allowable stress range across all component
                curves.
        """
        return np.maximum.reduce([curve.DS(N) for curve in self.curves])


class DoubleSlopeFatigueCurve(MultiSlopeFatigueCurve):
    r"""Bilinear (double-slope) Wöhler (S-N) fatigue curve.

    Implements a two-slope S-N curve commonly used in structural design standards
    (e.g., Eurocode 3, DNV, IEC standards). The curve transitions from slope $m_1$
    in the lower-cycle / higher-stress regime ($N \le N_{12}$) to slope $m_2$ in the
    higher-cycle / lower-stress regime ($N > N_{12}$) at an intersection knee point
    $(\Delta S_{12}, N_{12})$.

    Args:
        m1 (float): Logarithmic slope exponent for the lower-cycle regime ($N \le N_{12}$).
        m2 (float): Logarithmic slope exponent for the higher-cycle regime ($N > N_{12}$).
        DS12 (float): Stress range $\Delta S_{12}$ at the knee point where the two slopes meet.
        N12 (float): Number of cycles $N_{12}$ at the knee point where the two slopes meet.
    """

    def __init__(self, m1: float, m2: float, DS12: float, N12: float) -> None:
        r"""Initializes the DoubleSlopeFatigueCurve.

        Args:
            m1 (float): Logarithmic slope exponent for the lower-cycle regime ($N \le N_{12}$).
            m2 (float): Logarithmic slope exponent for the higher-cycle regime ($N > N_{12}$).
            DS12 (float): Stress range $\Delta S_{12}$ at the knee point where the two slopes meet.
            N12 (float): Number of cycles $N_{12}$ at the knee point where the two slopes meet.
        """
        curve1 = SingleSlopeFatigueCurve(m1, DS12, N12)
        curve2 = SingleSlopeFatigueCurve(m2, DS12, N12)
        super().__init__(curve1, curve2)


class BoltFatigueCurve(DoubleSlopeFatigueCurve):
    r"""Bolt fatigue S-N curve according to IEC 61400-6 AMD1.

    Constructs a design bilinear S-N curve for threaded bolts with logarithmic slopes
    $m_1 = 3$ (for $N \le 2 \times 10^6$) and $m_2 = 5$ (for $N > 2 \times 10^6$), with a
    knee point located at $N_{12} = 2 \times 10^6$ cycles.

    The characteristic reference stress range $\Delta \sigma_c$ at 2 million cycles
    is adjusted for bolt size effects according to IEC 61400-6 AMD1:

    - For $d \le 30\text{ mm}$:
      $$\Delta \sigma_c = \Delta S_{\text{ref}}$$
    - For $30\text{ mm} < d \le 72\text{ mm}$:
      $$\Delta \sigma_c = \Delta S_{\text{ref}} \cdot \left(\frac{0.030}{d}\right)^{0.1}$$
    - For $d > 72\text{ mm}$:
      $$\Delta \sigma_c = \Delta S_{\text{ref}} \cdot \left(\frac{0.030}{d}\right)^{0.1} \cdot \left(\frac{0.072}{d}\right)^{0.25}$$

    The knee point stress range for the design curve is obtained by dividing the
    characteristic strength by the material partial safety factor:

    $$\Delta S_{12} = \frac{\Delta \sigma_c}{\gamma_M}$$

    Args:
        diameter (float): Nominal bolt diameter $d$ in meters.
        DS_ref (float, optional): Reference characteristic stress range $\Delta S_{\text{ref}}$
            at $2 \times 10^6$ cycles before size and material factor corrections (in Pa).
            Defaults to 50 MPa (`50e6` Pa).
        gamma_M (float, optional): Material partial safety factor $\gamma_M$ for fatigue design.
            Defaults to 1.1.
    """

    def __init__(self, diameter: float, DS_ref: float = 50e6, gamma_M: float = 1.1) -> None:
        r"""Initializes the BoltFatigueCurve according to IEC 61400-6 AMD1.

        Args:
            diameter (float): Nominal bolt diameter $d$ in meters.
            DS_ref (float, optional): Reference characteristic stress range $\Delta S_{\text{ref}}$
                at $2 \times 10^6$ cycles (in Pa). Defaults to 50 MPa (`50e6` Pa).
            gamma_M (float, optional): Material partial safety factor $\gamma_M$.
                Defaults to 1.1.
        """
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
    r"""Fatigue analysis pipeline for a flange bolt.

    Evaluates the cyclic stress history, total cumulative fatigue damage, and expected
    fatigue life for a bolt in a flanged connection under external bending moment
    histories.

    The analysis transforms the flange shell moment Markov matrix into a bolt stress
    Markov matrix using the segment's non-linear transfer functions, accounting for
    axial force, bolt bending, bolt diameter effects, and dynamic stress multiplication.

    Attributes:
        fseg (FlangeSegment): The flange segment model containing geometry, bolt
            properties, preloads, and transfer functions.
        flange_mkvm (MarkovMatrix): The Markov matrix representing the external bending
            moment load history acting on the flange connection.
        custom_fatigue_curve (FatigueCurve | None): Optional custom S-N curve for the bolt.
            If `None`, an IEC 61400-6 AMD1 compliant `BoltFatigueCurve` is automatically
            constructed from the bolt nominal diameter. Defaults to `None`.
        allowable_damage (float): Maximum allowable cumulative fatigue damage limit
            $D_{\text{allow}}$. Defaults to 1.0.
        SMF (float): Stress Multiplication Factor applied to the bolt stress ranges to
            account for dynamic amplification or geometric uncertainties. Defaults to 1.0.
    """

    fseg: FlangeSegment
    flange_mkvm: MarkovMatrix
    custom_fatigue_curve: FatigueCurve | None = None
    allowable_damage: float = 1.0
    SMF: float = 1.0

    @functools.cached_property
    def fatigue_curve(self) -> FatigueCurve:
        """FatigueCurve: The S-N fatigue curve used for bolt damage evaluation."""
        return self.custom_fatigue_curve or BoltFatigueCurve(self.fseg.bolt.nominal_diameter)

    @functools.cached_property
    def bolt_mkvm(self) -> MarkovMatrix:
        """MarkovMatrix: The bolt stress Markov matrix including axial and bending stresses."""
        from math import log
        from .flangesegments import bolt_markov_matrix
        bending_factor = max(0.5, 0.5 + 0.5*log(self.fseg.bolt.nominal_diameter/0.036) / log(150/36))
        return bolt_markov_matrix(self.fseg, self.flange_mkvm, bending_factor, SMF=self.SMF)

    @functools.cached_property
    def damage(self) -> float:
        r"""float: Total cumulative fatigue damage $D = \sum \frac{n_i}{N(\Delta S_i)}$."""
        return self.fatigue_curve.cumulated_damage(self.bolt_mkvm)

    @functools.cached_property
    def fatigue_life(self) -> float:
        """float: Predicted fatigue life in the same time units as `flange_mkvm.duration`."""
        return self.allowable_damage / self.damage * self.flange_mkvm.duration
