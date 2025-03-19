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
The ``bolts`` module contains objects representing the bolt, washer and nut
fastener components. In particular it contains

- A ``MetricBolt`` class that generates generic bolts with metric screw thread and a
  ``StandardMetricBolt`` function that generates MetricBolt objects with standard properties.
- A ``FlatWasher`` class that generates generic flat washer and a ``ISOFlatWasher`` frunction
  that returns a ``FlatWasher`` with standard dimensions.
- A ``HexNut`` class that generates a generic hexagonal nut, a ``ISOHExNut`` function that
  generates a ``HexNut`` with ISO 4032 dimensions`` and a ``RoundNut`` function that generates
  a standard flanged ``HexNut``.
'''

from dataclasses import dataclass
from functools import cached_property, cache

from .logger import Logger, log_data
logger = Logger(__name__)

# UNITS OF MEASUREMENT
# Distance
m = 1
mm = 0.001*m
# Pressure
Pa = 1
MPa = 1e6*Pa
GPa = 1e9*Pa



class Bolt:
    pass


@dataclass
class BoltCrossSection:
    ''' Bolt circular cross-section.

    Args:
        diameter: cross-section diameter

    Attributes:
        diameter: cross-section diameter
    '''

    diameter: float

    @cached_property
    def area (self):
        ''' The cross-section area.'''
        from math import pi
        return pi * self.diameter**2 / 4

    @cached_property
    def second_moment_of_area (self):
        ''' The second moment of area of the cross-section.'''
        from math import pi
        return pi * self.diameter**4 / 64

    @cached_property
    def elastic_section_modulus (self):
        ''' The elastic section modulus of the cross-section.'''
        from math import pi
        return pi * self.diameter**3 / 32



@dataclass
class MetricBolt (Bolt):
    ''' Generic bolt with ISO 68-1 metric screw thread

    Args:
        nominal_diameter: The outermost diameter of the screw thread.

        thread_pitch: The pitch of the metric thread.

        yield_stress: Nominal yield stress (0.2% strain limit) of the bolt material.

        ultimate_tensile_stress: Nominal ultimate tensile stress of the bolt material.

        elastic_modulus: The Young's modulus of the bolt material.
            If omitted, it defaults to 210e9 N/m². Notice that the default value assumes
            that the chosen unit for distance is m and the chosen unit for forces is N. If
            that's not the case, you should enter the proper value of this parameter.

        poissons_ratio: The Poisson's ratio of the bolt material.
            If omitted, it defaults to 0.30.

        shank_length: The length of the shank. If omitted, it defaults to 0.

        shank_diameter_ratio: The ratio between the shank diameter and the bolt
            nominal diameter. If omitted, it defaults to 1, which means that the
            shank has the nominal diameter.

        stud: True if this is a stud bolt, False if it is not. If omitted, it defaults to False.

    The parameters must be expressed in a consistent system of units. For example,
    if you chose to input distances in mm and forces in N, then stresses must be
    expressed in N/mm². All the bolt attributes and methods will return values
    consistently with the input units of measurement.

    All the input parameters are also available as attributes of the generated
    object (e.g. ``bolt.shank_length``, ``bolt.yield_stress``, etc.).

    This instances of this calss are designed to be immutable, which means than
    changing an attribute after creating an object is not a good idea. If you
    need a different bolt with different attributes, create a new one.
    '''

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
    def designation (self):
        ''' Bolt designation string.

        For example, ``"M16"`` is the designation of a bolt with nominal
        diameter 16 mm.
        '''
        return f"M{int(self.nominal_diameter*1000)}"


    @cached_property
    def shank_diameter (self):
        ''' Diameter of the shank. '''
        return self.nominal_diameter * self.shank_diameter_ratio


    @cached_property
    def thread_height (self):
        ''' Height of the metric thread fundamental triangle (H).

        As defined in ISO 68-1:1998.
        '''
        return 0.5 * 3**0.5 * self.thread_pitch


    @cached_property
    def thread_basic_minor_diameter (self):
        ''' Basic minor diameter (d1).

        As defined in ISO 68-1:1998.'''
        return self.nominal_diameter - 2 * 5/8 * self.thread_height


    @cached_property
    def thread_basic_pitch_diameter (self):
        ''' Basic minor diameter (d2).

        As defined in ISO 68-1:1998.'''
        return self.nominal_diameter - 2 * 3/8 * self.thread_height


    @cached_property
    def thread_minor_diameter (self):
        ''' Minor diameter (d3).

        As defined in ISO 898-1:2013.'''
        return self.thread_basic_minor_diameter - self.thread_height/6


    @cached_property
    def nominal_cross_section (self):
        ''' Bolt cross-section with nominal diameter.

        Instance of `BoltCrossSection` class.
        '''
        return BoltCrossSection(self.nominal_diameter)


    @cached_property
    def shank_cross_section (self):
        ''' Bolt shank cross-section.

        Instance of `BoltCrossSection` class.
        '''
        return BoltCrossSection(self.shank_diameter)


    @cached_property
    def thread_cross_section (self):
        ''' Bolt cross-section used for tensile calculations.

        Instance of `BoltCrossSection` class.
        Ref. ISO 891-1:2013, section 9.1.6.1
        '''
        return BoltCrossSection(self.nominal_diameter - 13/12*self.thread_height)



    # --------------------------------------------------------------------------
    #   MATERIAL PROPERTIES
    # --------------------------------------------------------------------------

    @cached_property
    def shear_modulus (self):
        ''' Shear modulus G.

        Calculated from the Young's modulus and Poisson's ratio, under the
        assumption of isotropic and elastic bolt material.
        '''
        return 0.5 * self.elastic_modulus / (1 + self.poissons_ratio)



    # --------------------------------------------------------------------------
    #   MECHANICAL PROPERTIES
    # --------------------------------------------------------------------------

    def ultimate_tensile_capacity (self, standard="Eurocode"):
        ''' Evaluate the ultimate tensile force that the bolt can take.

        Args:
            standard (str): Standard according to which the ultimate tensile force
                should be calculated. Currently the only supported standard
                is *"Eurocode"* (EN 1993-1-8:2005).

        Returns:
            FRu (float): The bolt ultimate tensile force according to the specified
                standard.

        Raises:
            ValueError: if the requested standard is not supported.

        '''
        if standard.upper() == "EUROCODE":
            return 0.9 * self.ultimate_tensile_stress * self.thread_cross_section.area / 1.25
        else:
            raise ValueError(f"Unsupported standard: '{standard}'")


    def axial_stiffness (self, length):
        ''' Evaluate the axial stiffness of the bolt.

        Args:
            length (float): clamped length.

        Returns:
            Ka (float): axial stiffness of the bolt, according to VDI 2230,
                Part 1, Section 5.1.1.1.
        '''

        # Verify input validity
        assert length >= self.shank_length, "The lolt can't be shorter than its shank."

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


    def bending_stiffness (self, length):
        ''' Evaluates the bending stiffness of the bolt.

        Args:
            length (float): clamped length.

        Returns:
            Kb (float): bending stiffness of the bolt, according to VDI 2230,
                Part 1, Section 5.1.1.2.
        '''

        # Verify input validity
        assert length >= self.shank_length, "The lolt can't be shorter than its shank."

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
    def shank_cross_section_area (self):
        ''' Area of the shank transversal cross-section.

        **DEPRECATED**: use `bolt.shank_cross_section.area` instead.
        '''

        from .logger import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.shank_cross_section_area is deprecated; use MetricBolt.shank_cross_section.area instead.")

        from math import pi
        return pi * self.shank_diameter**2 / 4


    @cached_property
    def nominal_cross_section_area (self):
        ''' Area of a circle with nominal diameter.

        **DEPRECATED**: use `bolt.nominal_cross_section.area` instead
        '''

        from .logger import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.nominal_cross_section_area is deprecated; use MetricBolt.nominal_cross_section.area instead.")

        from math import pi
        return pi * self.nominal_diameter**2 / 4


    @cached_property
    def tensile_cross_section_area (self):
        ''' Tensile stress area, according to ISO 891-1:2013, section 9.1.6.1.

        **DEPRECATED**: use `bolt.thread_cross_section.area` instead
        '''

        from .logger import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.tensile_cross_section_area is deprecated; use MetricBolt.thread_cross_section.area instead.")

        from math import pi
        return pi * (self.nominal_diameter - 13/12*self.thread_height)**2 / 4


    @cached_property
    def tensile_moment_of_resistance (self):
        ''' Tensile moment of resistance, according to ISO 891-1:2013, section 9.1.6.1.

        **DEPRECATED**: use `bolt.thread_cross_section.elastic_section_modulus` instead
        '''

        from .logger import Logger
        logger = Logger(__name__)
        logger.warning("MetricBolt.tensile_moment_of_resistance is deprecated; use MetricBolt.thread_cross_section.elastic_section_modulus instead.")

        from math import pi
        return pi * (self.nominal_diameter - 13/12*self.thread_height) ** 3/32



def StandardMetricBolt (designation, material_grade, shank_length=0.0, shank_diameter_ratio=1.0, stud=False):
    ''' Create a metric bolt with standard dimensions.

    This function provides a convenient way for creating ``MetricBolt`` object,
    given the standard geometry designation (e.g. "M20") and the standard material
    grade designation (e.g. "8.8").

    Args:
        designation (str): The metric screw thread designation. The allowed values are:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20', 'M22',
            'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45', 'M48', 'M52', 'M56',
            'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

        material_grade (str): The material grade designation. The allowed values are:
            '4.6', '4.8', '5.6', '5.8', '6.8', '8.8', '9.8', '10.9' and '12.9' for
            carbon-steel bolts; 'A50', 'A70', 'A80' and 'A100' for austenitic bolts;
            'D70', 'D80' and 'D100' for duplex bolts; 'C50', 'C70', 'C80' and 'C110' for
            martensitic bolts; 'F45' and 'F60' for ferritic bolts.

        shank_length (float): The length of the shank.

        shank_diameter_ratio (float): The ratio between the shank diameter and the
            bolt nominal diameter.

        stud (bool): True if this is a stud bolt, False if it is not.

    Returns:
        bolt (MetricBolt): a MetricBolt instance with standard properties.
    '''

    geometry = _load_database('bolts.metric_screws')
    material = _load_database('bolts.materials')

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
    pass



@dataclass
class FlatWasher (Washer):
    ''' Generic flat washer.

    Args:
        outer_diameter: The outer diameter of the washer.

        inner_diameter: The hole diameter of the washer.

        elastic_modulus: The Young's modulus of the washer material.
            If omitted, it defaults to 210e9 N/m². Notice that the default value assumes
            that the chosen unit for distance is m and the chosen unit for forces is N. If
            that's not the case, you should enter the proper value of this parameter.

        poissons_ratio: The Poisson's ratio of the washer material.

    The parameters must be expressed in a consistent system of units. For example,
    if you chose to input distances in mm and forces in N, then stresses must be
    expressed in N/mm². All the bolt attributes and methods will return values
    consistently with the input units of measurement.

    All the input parameters are also available as attributes of the generated
    object (e.g. ``washer.thickness``, ``washer.poissons_ratio``, etc.).

    This instances of this calss are designed to be immutable, which means than
    changing an attribute after creating an object is not a good idea. If you
    need a different washer with different attributes, create a new one.
    '''

    outer_diameter: float
    inner_diameter: float
    thickness: float

    elastic_modulus: float = 210*GPa
    poissons_ratio: float = 0.3

    @cached_property
    def area (self):
        ''' Area of the washer flat surface '''
        from math import pi
        return pi/4 * (self.outer_diameter**2 - self.inner_diameter**2)

    @cached_property
    def axial_stiffness (self):
        ''' The compressive stiffness of the washer: t / EA'''
        return self.elastic_modulus * self.area / self.thickness



def ISOFlatWasher (designation):
    ''' Generates a standard washer according to ISO 7089.

    Args:
        designation (str): The metric screw thread designation. The allowed values are:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20', 'M22',
            'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45', 'M48', 'M52', 'M56',
            'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        washer (FlatWasher): a FlatWasher instance  having the standard dimensions
            defined in ISO 7089.

    For example, ``ISOFlatWasher("M16")`` will return a ``FlatWasher``
    instance with outer diameter 30 mm, hole diameter 17 mm and
    thickness 3 mm.
    '''

    params = _load_database("bolts.flat_washers")
    return FlatWasher(
        outer_diameter = params['outer_diameter'][designation],
        inner_diameter = params['hole_diameter'][designation],
        thickness = params['thickness'][designation])



class Nut:
    pass



@dataclass
class HexNut (Nut):
    ''' Generates a generic hexagonal nut.

    Args:
        nominal_diameter: The nominal diameter of the inner thread.

        thickness: The height of the bolt.

        inscribed_diameter: The diameter of the circle inscribed in the hexagon.
            Correponds to the distance between two opposite flats.

        circumscribed_diameter: The diameter of the circle circumscribed in the hexagon.
            Correponds to the distance between two opposite vertices..

        bearing_diameter: The outer diameter of the circular contatact surface
            between nut and washer.

        elastic_modulus: The Young's modulus of the nut material. If omitted, it
            defaults to 210e9 N/m². Notice that the default value assumes that the
            chosen unit for distance is m and the chosen unit for forces is N. If
            that's not the case, you should enter the proper value of this parameter.

        poissons_ratio: The Poisson's ratio of the nut material.

    The parameters must be expressed in a consistent system of units. For example,
    if you chose to input distances in mm and forces in N, then stresses must be
    expressed in N/mm². All the bolt attributes and methods will return values
    consistently with the input units of measurement.

    All the input parameters are also available as attributes of the generated
    object (e.g. ``washer.thickness``, ``washer.poissons_ratio``, etc.).

    This instances of this calss are designed to be immutable, which means than
    changing an attribute after creating an object is not a good idea. If you
    need a different nut with different attributes, create a new one.
    '''

    nominal_diameter: float         # nominal diameter of the thread
    thickness: float                # height of the nut
    inscribed_diameter: float       # distance between flats
    circumscribed_diameter: float   # distances between vertices

    bearing_diameter: float         # the diameter of the surface in contact with the washer

    elastic_modulus: float = 210*GPa
    poissons_ratio: float = 0.3


def ISOHexNut (designation):
    ''' Generates a standard Hex Nut.

    Args:
        designation (str): The metric screw thread designation. The allowed values are:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20', 'M22',
            'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45', 'M48', 'M52', 'M56',
            'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        nut (HexNut): a HexNut instance with dimensions according to ISO 4032.
    '''
    # params = _standard["hex_nuts"][designation]
    params = _load_database("bolts.hex_nuts")
    return HexNut(
        nominal_diameter = params["nominal_diameter"][designation],
        thickness = params["thickness"][designation],
        inscribed_diameter = params["inscribed_diameter"][designation],
        circumscribed_diameter = params["circumscribed_diameter"][designation],
        bearing_diameter = params["bearing_diameter"][designation]
    )


def RoundNut (designation):
    ''' Generates a standard round nut.


    Args:
        designation (str): The metric screw thread designation. The allowed values are:
            'M4', 'M5', 'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20', 'M22',
            'M24', 'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45', 'M48', 'M52', 'M56',
            'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

    Returns:
        nut (HexNut): a standard flanged nut.
    '''
    params = _standard["round_nuts"][designation]
    return HexNut(
        nominal_diameter = params["D_nom"],
        thickness = params["t"],
        inscribed_diameter = params["D_ins"],
        circumscribed_diameter = params["D_cir"],
        bearing_diameter = params["D_brg"]
    )


@cache
def _load_database (db_name):
    import os
    import pandas as pd
    db_path = os.path.join(os.path.dirname(__file__), f"data/{db_name}.csv")
    return pd.read_csv(db_path, index_col="id")



_standard = {


    "round_nuts": {

        "M36" : {"D_nom": 36*mm,     # nominal diameter
                 "t"    : 36*mm,     # nut thickness
                 "D_cir": 60.79*mm,     # nut circumscribed circle diameter
                 "D_ins": 55*mm,
                 "D_brg": 66*mm},    # nut inscribed circle diameter

        "M42" : {"D_nom":42*mm, "t":42.0*mm, "D_cir": 71.30*mm, "D_ins": 65*mm, "D_brg": 78*mm},
        "M48" : {"D_nom":48*mm, "t":48.0*mm, "D_cir": 82.60*mm, "D_ins": 75*mm, "D_brg": 92*mm},
        "M56" : {"D_nom":56*mm, "t":56.0*mm, "D_cir": 93.56*mm, "D_ins": 85*mm, "D_brg":105*mm},
        "M64" : {"D_nom":64*mm, "t":64.0*mm, "D_cir":104.86*mm, "D_ins": 95*mm, "D_brg":115*mm},
        "M72" : {"D_nom":72*mm, "t":72.0*mm, "D_cir":116.16*mm, "D_ins":105*mm, "D_brg":125*mm},
        "M80" : {"D_nom":80*mm, "t":80.0*mm, "D_cir":127.46*mm, "D_ins":115*mm, "D_brg":140*mm}
    }
}
