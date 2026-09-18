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

"""Bolt, washer, and nut fastener components.

This module contains classes and factory functions representing bolt, washer,
and nut fastener components. In particular it provides:

- A `MetricBolt` class that models generic bolts with metric screw threads and a
  `StandardMetricBolt` factory function that generates `MetricBolt` objects with
  standard properties.
- A `FlatWasher` class that models generic flat washers and an `ISOFlatWasher`
  factory function that returns a `FlatWasher` with ISO 7089 standard dimensions.
- A `HexNut` class that models generic hexagonal nuts, an `ISOHexNut` factory
  function that generates a `HexNut` with ISO 4032 dimensions, and a `RoundNut`
  factory function that generates a standard flanged nut.
"""

from dataclasses import dataclass
from functools import cached_property

from .utils import Logger, log_data
logger = Logger(__name__)

from .utils import load_csv_database


# UNITS OF MEASUREMENT
# Distance
m = 1
mm = 0.001*m
# Pressure
Pa = 1
MPa = 1e6*Pa
GPa = 1e9*Pa



class Bolt:
    """Base class for bolt representations."""
    pass


@dataclass
class BoltCrossSection:
    """Circular cross-section of a bolt.

    Attributes:
        diameter (float): Cross-section diameter.
    """

    diameter: float

    @cached_property
    def area(self):
        """float: Cross-sectional area."""
        from math import pi
        return pi * self.diameter**2 / 4

    @cached_property
    def second_moment_of_area(self):
        """float: Second moment of area of the circular cross-section."""
        from math import pi
        return pi * self.diameter**4 / 64

    @cached_property
    def elastic_section_modulus(self):
        """float: Elastic section modulus of the circular cross-section."""
        from math import pi
        return pi * self.diameter**3 / 32



@dataclass
class MetricBolt(Bolt):
    """Generic bolt with ISO 68-1 metric screw thread.

    The parameters must be expressed in a consistent system of units. For
    example, if you choose to input distances in meters (m) and forces in
    newtons (N), then stresses must be expressed in pascals (N/m²). All bolt
    attributes and methods return values consistent with the input units of
    measurement.

    All input parameters are also available as attributes of the generated
    instance (e.g. `bolt.shank_length`, `bolt.yield_stress`, etc.).

    Instances of this class are designed to be immutable; changing attributes
    after instantiation is not recommended. If a bolt with different parameters
    is needed, instantiate a new one.

    Attributes:
        nominal_diameter (float): Outermost (nominal) diameter of the screw thread.
        thread_pitch (float): Pitch of the metric thread.
        yield_stress (float): Nominal yield stress (0.2% strain limit) of the bolt material.
        ultimate_tensile_stress (float): Nominal ultimate tensile stress of the bolt material.
        elastic_modulus (float): Young's modulus of the bolt material.
            Defaults to 210 GPa (210e9 N/m²).
        poissons_ratio (float): Poisson's ratio of the bolt material.
            Defaults to 0.3.
        shank_length (float): Length of the unthreaded shank. Defaults to 0.0.
        shank_diameter_ratio (float): Ratio between the shank diameter and the
            nominal diameter. Defaults to 1.0 (shank has nominal diameter).
        stud (bool): True if this is a stud bolt, False otherwise. Defaults to False.
    """

    nominal_diameter: float
    thread_pitch: float

    yield_stress: float
    ultimate_tensile_stress: float
    elastic_modulus: float = 210*GPa
    poissons_ratio: float = 0.3

    shank_length: float = 0
    shank_diameter_ratio: float = 1
    stud: bool = False


    # --------------------------------------------------------------------------
    #   GEOMETRY
    # --------------------------------------------------------------------------

    @cached_property
    def designation(self):
        """str: Bolt designation string (e.g. 'M16' for a bolt with 16 mm diameter)."""
        return f"M{int(self.nominal_diameter*1000)}"


    @cached_property
    def shank_diameter(self):
        """float: Diameter of the unthreaded shank."""
        return self.nominal_diameter * self.shank_diameter_ratio


    @cached_property
    def thread_height(self):
        """float: Height of the metric thread fundamental triangle (H), per ISO 68-1:1998."""
        return 0.5 * 3**0.5 * self.thread_pitch


    @cached_property
    def thread_basic_minor_diameter(self):
        """float: Basic minor diameter (d1), per ISO 68-1:1998."""
        return self.nominal_diameter - 2 * 5/8 * self.thread_height


    @cached_property
    def thread_basic_pitch_diameter(self):
        """float: Basic pitch diameter (d2), per ISO 68-1:1998."""
        return self.nominal_diameter - 2 * 3/8 * self.thread_height


    @cached_property
    def thread_minor_diameter(self):
        """float: Minor diameter (d3), per ISO 898-1:2013."""
        return self.thread_basic_minor_diameter - self.thread_height/6


    @cached_property
    def nominal_cross_section(self):
        """BoltCrossSection: Bolt cross-section with nominal diameter."""
        return BoltCrossSection(self.nominal_diameter)


    @cached_property
    def shank_cross_section(self):
        """BoltCrossSection: Bolt shank cross-section."""
        return BoltCrossSection(self.shank_diameter)


    @cached_property
    def thread_cross_section(self):
        """BoltCrossSection: Bolt cross-section used for tensile calculations (ISO 898-1:2013, sec. 9.1.6.1)."""
        return BoltCrossSection(self.nominal_diameter - 13/12*self.thread_height)



    # --------------------------------------------------------------------------
    #   MATERIAL PROPERTIES
    # --------------------------------------------------------------------------

    @cached_property
    def shear_modulus(self):
        """float: Shear modulus G, under the assumption of isotropic linear elastic material."""
        return 0.5 * self.elastic_modulus / (1 + self.poissons_ratio)



    # --------------------------------------------------------------------------
    #   MECHANICAL PROPERTIES
    # --------------------------------------------------------------------------

    def ultimate_tensile_capacity(self, standard="Eurocode"):
        """Evaluate the ultimate tensile force that the bolt can sustain.

        Args:
            standard (str): Standard according to which the ultimate tensile force
                should be calculated. Currently supported: "Eurocode" (EN 1993-1-8:2005).
                Defaults to "Eurocode".

        Returns:
            float: The bolt ultimate tensile force according to the specified standard.

        Raises:
            ValueError: If the requested standard is not supported.
        """
        if standard.upper() == "EUROCODE":
            return 0.9 * self.ultimate_tensile_stress * self.thread_cross_section.area / 1.25
        else:
            raise ValueError(f"Unsupported standard: '{standard}'")


    def axial_stiffness(self, length):
        """Evaluate the axial stiffness of the clamped bolt.

        Calculates the axial stiffness according to VDI 2230 Part 1,
        Section 5.1.1.1.

        Args:
            length (float): Clamped length.

        Returns:
            float: Axial stiffness of the bolt.

        Raises:
            AssertionError: If `length` is less than `shank_length`.
        """

        # Verify input validity
        assert length >= self.shank_length, "The bolt can't be shorter than its shank."

        # Common variables
        from math import pi
        E = self.elastic_modulus
        An = self.nominal_cross_section.area
        As = self.shank_cross_section.area
        At = pi * self.thread_minor_diameter**2 / 4

        # Resilience of unthreaded part
        L1 = self.shank_length
        d1 = L1 / (E * As)

        # Resilience at the minor diameter of the engaged bolt thread
        LG = 0.5 * self.nominal_diameter
        dG = LG / (E * At)

        # Resilience of the nut
        LM = 0.4 * self.nominal_diameter
        dM = LM / (E * An)

        # Resilience of threaded part
        LGew = length - self.shank_length
        dGew = LGew / (E * At)

        # Resilience of hex head
        LSK = 0.5 * self.nominal_diameter
        dSK = LSK / (E * An)

        # Total stiffness
        if self.stud:
            return 1 / (d1 + 2*dG + 2*dM + dGew)
        else:
            return 1 / (d1 + dG + dM + dGew + dSK)


    def bending_stiffness(self, length):
        """Evaluate the bending stiffness of the clamped bolt.

        Calculates the bending stiffness according to VDI 2230 Part 1,
        Section 5.1.1.2.

        Args:
            length (float): Clamped length.

        Returns:
            float: Bending stiffness of the bolt.

        Raises:
            AssertionError: If `length` is less than `shank_length`.
        """

        # Verify input validity
        assert length >= self.shank_length, "The bolt can't be shorter than its shank."

        # Common variables
        from math import pi
        E = self.elastic_modulus
        In = pi * self.nominal_diameter**4 / 64
        Is = pi * self.shank_diameter**4 / 64
        It = pi * self.thread_minor_diameter**4 / 64

        # Bending resilience of unthreaded part
        L1 = self.shank_length
        b1 = L1 / (E * Is)

        # Bending resilience at the minor diameter of the engaged bolt thread
        LG = 0.5 * self.nominal_diameter
        bG = LG / (E * It)

        # Bending resilience of the nut
        LM = 0.4 * self.nominal_diameter
        bM = LM / (E * In)

        # Bending resilience of threaded part
        LGew = length - self.shank_length
        bGew = LGew / (E * It)

        # Bending resilience of hex head
        LSK = 0.5 * self.nominal_diameter
        bSK = LSK / (E * In)

        log_data(self, beta_sk=bSK, beta_1=b1, beta_Gew=bGew, beta_G=bG, beta_M=bM)

        # Total bending stiffness
        if self.stud:
            return 1 / (b1 + 2*bG + 2*bM + bGew)
        else:
            return 1 / (b1 + bG + bM + bGew + bSK)



    # --------------------------------------------------------------------------
    #   DEPRECATED ATTRIBUTES AND METHODS
    # --------------------------------------------------------------------------

    @cached_property
    def shank_cross_section_area(self):
        """float: Area of the shank transversal cross-section.

        Deprecated:
            Use `MetricBolt.shank_cross_section.area` instead.
        """

        from .utils import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.shank_cross_section_area is deprecated; use MetricBolt.shank_cross_section.area instead.")

        from math import pi
        return pi * self.shank_diameter**2 / 4


    @cached_property
    def nominal_cross_section_area(self):
        """float: Area of a circle with nominal diameter.

        Deprecated:
            Use `MetricBolt.nominal_cross_section.area` instead.
        """

        from .utils import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.nominal_cross_section_area is deprecated; use MetricBolt.nominal_cross_section.area instead.")

        from math import pi
        return pi * self.nominal_diameter**2 / 4


    @cached_property
    def tensile_cross_section_area(self):
        """float: Tensile stress area, according to ISO 898-1:2013, section 9.1.6.1.

        Deprecated:
            Use `MetricBolt.thread_cross_section.area` instead.
        """

        from .utils import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.tensile_cross_section_area is deprecated; use MetricBolt.thread_cross_section.area instead.")

        from math import pi
        return pi * (self.nominal_diameter - 13/12*self.thread_height)**2 / 4


    @cached_property
    def tensile_moment_of_resistance(self):
        """float: Tensile moment of resistance, according to ISO 898-1:2013, section 9.1.6.1.

        Deprecated:
            Use `MetricBolt.thread_cross_section.elastic_section_modulus` instead.
        """

        from .utils import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.tensile_moment_of_resistance is deprecated; use MetricBolt.thread_cross_section.elastic_section_modulus instead.")

        from math import pi
        return pi * (self.nominal_diameter - 13/12*self.thread_height) ** 3/32



def StandardMetricBolt(designation, material_grade, shank_length=0.0, shank_diameter_ratio=1.0, stud=False):
    """Create a metric bolt with standard dimensions and material properties.

    This function provides a convenient way of creating a `MetricBolt` object,
    given the standard geometry designation (e.g. "M20") and the standard material
    grade designation (e.g. "8.8").

    Args:
        designation (str): Metric screw thread designation. Allowed values:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20',
            'M22', 'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45',
            'M48', 'M52', 'M56', 'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.
        material_grade (str): Material grade designation. Allowed values:
            - Carbon-steel: '4.6', '4.8', '5.6', '5.8', '6.8', '8.8', '9.8',
              '10.9', '12.9'
            - Austenitic stainless: 'A50', 'A70', 'A80', 'A100'
            - Duplex stainless: 'D70', 'D80', 'D100'
            - Martensitic stainless: 'C50', 'C70', 'C80', 'C110'
            - Ferritic stainless: 'F45', 'F60'
        shank_length (float, optional): Length of the unthreaded shank.
            Defaults to 0.0.
        shank_diameter_ratio (float, optional): Ratio between the shank diameter
            and the bolt nominal diameter. Defaults to 1.0.
        stud (bool, optional): True if this is a stud bolt, False otherwise.
            Defaults to False.

    Returns:
        MetricBolt: A `MetricBolt` instance with standard properties.
    """

    geometry = load_csv_database('bolts.metric_screws')
    material = load_csv_database('bolts.materials')

    return MetricBolt(
        nominal_diameter = geometry['nominal_diameter'][designation],
        thread_pitch = geometry['course_pitch'][designation],
        yield_stress = material['yield_stress'][material_grade],
        ultimate_tensile_stress = material['ultimate_tensile_stress'][material_grade],
        elastic_modulus = material['youngs_modulus'][material_grade],
        poissons_ratio = material['poissons_ratio'][material_grade],
        shank_length = shank_length,
        shank_diameter_ratio = shank_diameter_ratio,
        stud = stud)



class Washer:
    """Base class for washer representations."""
    pass



@dataclass
class FlatWasher(Washer):
    """Generic flat washer.

    The parameters must be expressed in a consistent system of units. For
    example, if you choose to input distances in meters (m) and forces in
    newtons (N), then stresses must be expressed in pascals (N/m²). All
    attributes and methods return values consistent with the input units of
    measurement.

    All input parameters are also available as attributes of the generated
    instance (e.g. `washer.thickness`, `washer.poissons_ratio`, etc.).

    Instances of this class are designed to be immutable; changing attributes
    after instantiation is not recommended. If a washer with different
    attributes is needed, instantiate a new one.

    Attributes:
        outer_diameter (float): Outer diameter of the washer.
        inner_diameter (float): Inner (hole) diameter of the washer.
        thickness (float): Thickness of the washer.
        elastic_modulus (float): Young's modulus of the washer material.
            Defaults to 210 GPa (210e9 N/m²).
        poissons_ratio (float): Poisson's ratio of the washer material.
            Defaults to 0.3.
    """

    outer_diameter: float
    inner_diameter: float
    thickness: float

    elastic_modulus: float = 210*GPa
    poissons_ratio: float = 0.3

    @cached_property
    def area(self):
        """float: Surface area of the flat annular face."""
        from math import pi
        return pi/4 * (self.outer_diameter**2 - self.inner_diameter**2)

    @cached_property
    def axial_stiffness(self):
        """float: Compressive axial stiffness of the washer (EA / t)."""
        return self.elastic_modulus * self.area / self.thickness



def ISOFlatWasher(designation):
    """Generate a standard flat washer according to ISO 7089.

    Args:
        designation (str): Metric screw thread designation. Allowed values:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20',
            'M22', 'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45',
            'M48', 'M52', 'M56', 'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        FlatWasher: A `FlatWasher` instance having standard dimensions defined
            in ISO 7089 (e.g., for "M16", outer diameter 30 mm, hole diameter
            17 mm, thickness 3 mm).
    """

    params = load_csv_database("bolts.flat_washers")
    return FlatWasher(
        outer_diameter = params['outer_diameter'][designation],
        inner_diameter = params['hole_diameter'][designation],
        thickness = params['thickness'][designation])



class Nut:
    """Base class for nut representations."""
    pass



@dataclass
class HexNut(Nut):
    """Generic hexagonal nut.

    The parameters must be expressed in a consistent system of units. For
    example, if you choose to input distances in meters (m) and forces in
    newtons (N), then stresses must be expressed in pascals (N/m²). All
    attributes and methods return values consistent with the input units of
    measurement.

    All input parameters are also available as attributes of the generated
    instance (e.g. `nut.thickness`, `nut.bearing_diameter`, etc.).

    Instances of this class are designed to be immutable; changing attributes
    after instantiation is not recommended. If a nut with different attributes
    is needed, instantiate a new one.

    Attributes:
        nominal_diameter (float): Nominal diameter of the inner thread.
        thickness (float): Height/thickness of the nut.
        inscribed_diameter (float): Diameter of the circle inscribed in the
            hexagon (distance between opposite flats).
        circumscribed_diameter (float): Diameter of the circle circumscribed
            around the hexagon (distance between opposite vertices).
        bearing_diameter (float): Outer diameter of the circular contact surface
            between nut and washer or flange.
        elastic_modulus (float): Young's modulus of the nut material.
            Defaults to 210 GPa (210e9 N/m²).
        poissons_ratio (float): Poisson's ratio of the nut material.
            Defaults to 0.3.
    """

    nominal_diameter: float         # nominal diameter of the thread
    thickness: float                # height of the nut
    inscribed_diameter: float       # distance between flats
    circumscribed_diameter: float   # distances between vertices

    bearing_diameter: float         # the diameter of the surface in contact with the washer

    elastic_modulus: float = 210*GPa
    poissons_ratio: float = 0.3


def ISOHexNut(designation):
    """Generate a standard hexagonal nut according to ISO 4032.

    Args:
        designation (str): Metric screw thread designation. Allowed values:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20',
            'M22', 'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45',
            'M48', 'M52', 'M56', 'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        HexNut: A `HexNut` instance with dimensions according to ISO 4032.
    """
    params = load_csv_database("bolts.hex_nuts")
    return HexNut(
        nominal_diameter = params["nominal_diameter"][designation],
        thickness = params["thickness"][designation],
        inscribed_diameter = params["inscribed_diameter"][designation],
        circumscribed_diameter = params["circumscribed_diameter"][designation],
        bearing_diameter = params["bearing_diameter"][designation]
    )


def RoundNut(designation):
    """Generate a standard flanged round nut.

    Args:
        designation (str): Metric screw thread designation. Allowed values:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20',
            'M22', 'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45',
            'M48', 'M52', 'M56', 'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        HexNut: A `HexNut` instance configured with dimensions of a standard
            flanged round nut.
    """
    params = load_csv_database("bolts.round_nuts")
    return HexNut(
        nominal_diameter = params["nominal_diameter"][designation],
        thickness = params["thickness"][designation],
        inscribed_diameter = params["inscribed_diameter"][designation],
        circumscribed_diameter = params["circumscribed_diameter"][designation],
        bearing_diameter = params["bearing_diameter"][designation]
    )
