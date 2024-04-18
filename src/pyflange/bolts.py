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
The ``bolts`` module contains a ``MetricBolt`` class that generates generic
bolts with metric screw thread and a ``StandardMetricBolt`` function that
generates MetricBolt objects with standard properties.
'''

from dataclasses import dataclass
from functools import cached_property

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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
class MetricBolt (Bolt):
    ''' Generates a generic bolt with metric screw thread accoridng ISO 68-1, defined
    by the following parameters:

    - ``nominal_diameter`` : ``float``
        The outermost diameter of the screw thread.

    - ``thread_pitch`` : ``float``
        The pitch of the metric thread.

    - ``yield_stress`` : ``float``
        Nominal yield stress (0.2% strain limit) of the bolt material.

    - ``ultimate_tensile_stress`` : ``float``
        Nominal ultimate tensile stress of the bolt material.

    - ``elastic_modulus`` : ``float`` [optional]
        The Young's modulus of the bolt material.
        If omitted, it defaults to 210e9 N/m². Notice that the default value assumes
        that the chosen unit for distance is m and the chosen unit for forces is N. If
        that's not the case, you should enter the proper value of this parameter.

    - ``poissons_ratio`` : ``float`` [optional]
        The Poisson's ratio of the bolt material.
        If omitted, it defaults to 0.30.

    - ``shank_length`` : ``float`` [optional]
        The length of the shank. If omitted, it defaults to 0.

    - ``shank_diameter_ratio`` : ``float`` [optional]
        The ratio between the shank diameter and the bolt nominal diameter.
        If omitted, it defaults to 1, which means that the shank hs the
        nominal diameter.

    - ``stud`` : ``bool`` [optional]
        True if this is a stud bolt, False if it is not.
        If omitted, it defaults to False.

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


    # -------------------------------------------------------------
    #   GEOMETRY
    # -------------------------------------------------------------

    @cached_property
    def designation (self):
        ''' Bolt designation string, which is, for example, ``"M16"``
        for a bolt with nominal diameter 16 mm.
        '''
        return f"M{int(self.nominal_diameter*1000)}"


    @cached_property
    def shank_diameter (self):
        ''' Diameter of the shank. '''
        return self.nominal_diameter * self.shank_diameter_ratio


    @cached_property
    def thread_height (self):
        ''' Height of the metric thread fundamental triangle (H),
        as defined in ISO 68-1:1998.
        '''
        return 0.5 * 3**0.5 * self.thread_pitch


    @cached_property
    def thread_basic_minor_diameter (self):
        ''' Basic minor diameter (d1) as defined in ISO 68-1:1998.'''
        return self.nominal_diameter - 2 * 5/8 * self.thread_height


    @cached_property
    def thread_basic_pitch_diameter (self):
        ''' Basic minor diameter (d2) as defined in ISO 68-1:1998.'''
        return self.nominal_diameter - 2 * 3/8 * self.thread_height


    @cached_property
    def thread_minor_diameter (self):
        ''' Minor diameter (d3) as defined in ISO 898-1:2013.'''
        return self.thread_basic_minor_diameter - self.thread_height/6


    @cached_property
    def shank_cross_section_area (self):
        ''' Area of the shank transversal cross-section.'''
        from math import pi
        return pi * self.shank_diameter**2 / 4


    @cached_property
    def nominal_cross_section_area (self):
        ''' Area of a circle with nominal diameter.'''
        from math import pi
        return pi * self.nominal_diameter**2 / 4


    @cached_property
    def tensile_cross_section_area (self):
        ''' Tensile stress area.'''
        from math import pi
        return pi * (self.nominal_diameter - 1.081*self.thread_height)**2 / 4



    # -------------------------------------------------------------
    #   MATERIAL PROPERTIES
    # -------------------------------------------------------------

    @cached_property
    def shear_modulus (self):
        ''' Shear modulus G, calculated from the Young's modulus and
        Poisson's ratio, under the assumption of isotropic and elastic
        bolt material.'''
        return 0.5 * self.elastic_modulus / (1 + self.poissons_ratio)



    # -------------------------------------------------------------
    #   MECHANICAL PROPERTIES
    # -------------------------------------------------------------

    def ultimate_tensile_capacity (self, standard="Eurocode"):
        ''' Returns the design ultimate tensile force that the bolt can take,
        according to a given standard.

        Currently the only standard available is "Eurocode", which is also the
        default value of the ``standard`` parameter.
        '''
        if standard == "Eurocode":
            return 0.9 * self.ultimate_tensile_stress * self.tensile_cross_section_area / 1.25
        else:
            raise ValueError(f"Unsupported standard: '{standard}'")


    def axial_stiffness (self, length):
        ''' Given a clamped ``length``, returns the axial stiffness of the bolt,
        according to VDI 2230, Part 1, Section 5.1.1.1.
        '''

        # Verify input validity
        assert length >= self.shank_length, "The lolt can't be shorter than its shank."

        # Common variables
        from math import pi
        E = self.elastic_modulus
        An = self.nominal_cross_section_area
        As = self.shank_cross_section_area
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
        ''' Given a clamped ``length``, returns the bending stiffness of the bolt,
        according to VDI 2230, Part 1, Section 5.1.1.2.
        '''

        logger.debug(f"{self.designation} BENDING STIFFNESS FOR BOLT LENGTH L = {length*1000:.1f} mm")

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

        logger.debug(f"Bending resilience of hex head: beta_Sk = {bSK*1e9} rad/(GN.m)")
        logger.debug(f"Bending resilience of unthreaded part beta_1 = {b1*1e9} rad/(GN.m)")
        logger.debug(f"Bending resilience of threaded part: beta_Gew = {bGew*1e9} rad/(GN.m)")
        logger.debug(f"Bending resilience of minor diameter of engaged bolt thread: beta_G = {bG*1e9} rad/(GN.m)")
        logger.debug(f"Bending resilience of nut: beta_M = {bM*1e9} rad/(GN.m)")

        # Total bending stiffness
        if self.stud:
            return 1 / (b1 + 2*bG + 2*bM + bGew)
        else:
            return 1 / (b1 + bG + bM + bGew + bSK)




def StandardMetricBolt (designation, material_grade, shank_length=0.0, shank_diameter_ratio=1.0, stud=False):
    ''' This function provides a convenient way for creating ``MetricBolt``
    object, given the standard geometry designation (e.g. "M20") and the
    standard material grade designation (e.g. "8.8").

    The required parameters are:

    - ``designation`` : ``str``
        The metric screw thread designation. The allowed values are: 'M4', 'M5',
        'M6', 'M8', 'M10', 'M12', 'M14', 'M16', 'M18', 'M20', 'M22', 'M24',
        'M27', 'M30', 'M33', 'M36', 'M39', 'M42', 'M45', 'M48', 'M52', 'M56',
        'M60', 'M64', 'M72', 'M80', 'M90', 'M100'.

        This parameter corresponds to a standard nominal diameter and a relevant
        thread pitch value (i.e. the standard coarse pitch value).

    - ``material_grade`` : ``str``
        The material grade designation. The allowed values are: '4.6', '4.8',
        '5.6', '5.8', '6.8', '8.8', '9.8', '10.9' and '12.9' for carbon-steel
        bolts; 'A50', 'A70', 'A80' and 'A100' for austenitic bolts; 'D70', 'D80'
        and 'D100' for duplex bolts; 'C50', 'C70', 'C80' and 'C110' for
        martensitic bolts; 'F45' and 'F60' for ferritic bolts.

        This parameter corresponds to a standard set of the paraters ``yield_stress``,
        ``ultimate_tensile_stress``, ``elastic_modulus`` and ``poissons_ratio``.

    - ``shank_length`` : ``float`` [optional]
        The length of the shank. If omitted, it defaults to 0.

    - ``shank_diameter_ratio`` : ``float`` [optional]
        The ratio between the shank diameter and the bolt nominal diameter.
        If omitted, it defaults to 1, which means that the shank hs the
        nominal diameter.

    - ``stud`` : ``bool`` [optional]
        True if this is a stud bolt, False if it is not.
        If omitted, it defaults to False.

    '''

    geometry = _standard['geometry'][designation]
    meterial = _standard['materials'][material_grade]

    return MetricBolt(
        nominal_diameter = geometry['D'],
        thread_pitch = geometry['Pc'],
        yield_stress = material['fy'],
        ultimate_tensile_stress = material['fu'],
        elastic_modulus = material['E'],
        poissons_ratio = material['nu'],
        shank_length = shank_length,
        shank_diameter_ratio = shank_diameter_ratio,
        stud = stud)

_standard = {

    "geometry": {

        "M4"  : {"D"    : 4.00*mm,     # nominal diameter
                "Pc"   : 0.70*mm,     # coarse screw thead pitch
                "Pf"   : None,        # fine screaw thead pitch
                "t_nut": 3.20*mm,     # nut thickness
                "t_hex": 2.80*mm,     # hex head thickness
                "D_hex": 7.66*mm,     # hex head/nut circumscribed circle diameter
                "d_hex": 7.00*mm,     # hex head/nut inscribed circle diameter
                "d_was": 4.30*mm,     # washer hole diameter ISO 7089
                "D_was": 9.00*mm,     # washer outer diameter ISO 7089
                "t_was": 0.80*mm},    # washer thickness ISO 7089

        "M5"  : {"D": 5*mm, "Pc":0.80*mm, "Pf":None, "t_nut": 4.7*mm, "t_hex": 3.5*mm, "D_hex":  8.79*mm, "d_hex":  8*mm, "d_was":  5.3*mm, "D_was": 10*mm, "t_was": 1.0*mm},
        "M6"  : {"D": 6*mm, "Pc":1.00*mm, "Pf":None, "t_nut": 5.2*mm, "t_hex": 4.0*mm, "D_hex": 11.05*mm, "d_hex": 10*mm, "d_was":  6.4*mm, "D_was": 12*mm, "t_was": 1.6*mm},
        "M8"  : {"D": 8*mm, "Pc":1.25*mm, "Pf":1.00, "t_nut": 6.8*mm, "t_hex": 5.3*mm, "D_hex": 14.38*mm, "d_hex": 13*mm, "d_was":  8.4*mm, "D_was": 16*mm, "t_was": 1.6*mm},
        "M10" : {"D":10*mm, "Pc":1.50*mm, "Pf":1.25, "t_nut": 8.4*mm, "t_hex": 6.4*mm, "D_hex": 18.90*mm, "d_hex": 17*mm, "d_was": 10.5*mm, "D_was": 20*mm, "t_was": 2.0*mm},
        "M12" : {"D":12*mm, "Pc":1.75*mm, "Pf":1.50, "t_nut":10.8*mm, "t_hex": 7.5*mm, "D_hex": 21.10*mm, "d_hex": 19*mm, "d_was": 13.0*mm, "D_was": 24*mm, "t_was": 2.5*mm},
        "M14" : {"D":14*mm, "Pc":2.00*mm, "Pf":1.50, "t_nut":12.8*mm, "t_hex": 8.8*mm, "D_hex": 24.49*mm, "d_hex": 22*mm, "d_was": 15.0*mm, "D_was": 28*mm, "t_was": 2.5*mm},  # 2nd choice
        "M16" : {"D":16*mm, "Pc":2.00*mm, "Pf":1.50, "t_nut":14.8*mm, "t_hex":10.0*mm, "D_hex": 26.75*mm, "d_hex": 24*mm, "d_was": 17.0*mm, "D_was": 30*mm, "t_was": 3.0*mm},
        "M18" : {"D":18*mm, "Pc":2.50*mm, "Pf":2.00, "t_nut":15.8*mm, "t_hex":11.5*mm, "D_hex": 30.14*mm, "d_hex": 27*mm, "d_was": 19.0*mm, "D_was": 34*mm, "t_was": 3.0*mm},  # 2nd choice
        "M20" : {"D":20*mm, "Pc":2.50*mm, "Pf":2.00, "t_nut":18.0*mm, "t_hex":12.5*mm, "D_hex": 33.53*mm, "d_hex": 30*mm, "d_was": 21.0*mm, "D_was": 37*mm, "t_was": 3.0*mm},
        "M22" : {"D":22*mm, "Pc":2.50*mm, "Pf":2.00, "t_nut":19.4*mm, "t_hex":14.0*mm, "D_hex": 35.72*mm, "d_hex": 32*mm, "d_was": 23.0*mm, "D_was": 39*mm, "t_was": 3.0*mm},  # 2nd choice
        "M24" : {"D":24*mm, "Pc":3.00*mm, "Pf":2.00, "t_nut":21.5*mm, "t_hex":15.0*mm, "D_hex": 39.98*mm, "d_hex": 36*mm, "d_was": 25.0*mm, "D_was": 44*mm, "t_was": 4.0*mm},
        "M27" : {"D":27*mm, "Pc":3.00*mm, "Pf":2.00, "t_nut":23.8*mm, "t_hex":17.0*mm, "D_hex": 45.20*mm, "d_hex": 41*mm, "d_was": 28.0*mm, "D_was": 50*mm, "t_was": 4.0*mm},  # 2nd choice
        "M30" : {"D":30*mm, "Pc":3.50*mm, "Pf":2.00, "t_nut":25.6*mm, "t_hex":18.7*mm, "D_hex": 50.85*mm, "d_hex": 46*mm, "d_was": 31.0*mm, "D_was": 56*mm, "t_was": 4.0*mm},
        "M33" : {"D":33*mm, "Pc":3.50*mm, "Pf":2.00, "t_nut":28.7*mm, "t_hex":21.0*mm, "D_hex": 55.37*mm, "d_hex": 50*mm, "d_was": 34.0*mm, "D_was": 60*mm, "t_was": 5.0*mm},  # 2nd choice
        "M36" : {"D":36*mm, "Pc":4.00*mm, "Pf":3.00, "t_nut":31.0*mm, "t_hex":22.5*mm, "D_hex": 60.79*mm, "d_hex": 55*mm, "d_was": 37.0*mm, "D_was": 66*mm, "t_was": 5.0*mm},
        "M39" : {"D":39*mm, "Pc":4.00*mm, "Pf":3.00, "t_nut":33.4*mm, "t_hex":25.0*mm, "D_hex": 66.44*mm, "d_hex": 60*mm, "d_was": 40.0*mm, "D_was": 72*mm, "t_was": 6.0*mm},  # 2nd choice
        "M42" : {"D":42*mm, "Pc":4.50*mm, "Pf":3.00, "t_nut":34.0*mm, "t_hex":26.0*mm, "D_hex": 71.30*mm, "d_hex": 65*mm, "d_was": 43.0*mm, "D_was": 78*mm, "t_was": 7.0*mm},
        "M45" : {"D":45*mm, "Pc":4.50*mm, "Pf":3.00, "t_nut":36.0*mm, "t_hex":28.0*mm, "D_hex": 76.95*mm, "d_hex": 70*mm, "d_was": 46.0*mm, "D_was": 85*mm, "t_was": 7.0*mm},  # 2nd choice
        "M48" : {"D":48*mm, "Pc":5.00*mm, "Pf":3.00, "t_nut":38.0*mm, "t_hex":30.0*mm, "D_hex": 82.60*mm, "d_hex": 75*mm, "d_was": 50.0*mm, "D_was": 92*mm, "t_was": 8.0*mm},
        "M52" : {"D":52*mm, "Pc":5.00*mm, "Pf":4.00, "t_nut":42.0*mm, "t_hex":33.0*mm, "D_hex": 88.25*mm, "d_hex": 80*mm, "d_was": 54.0*mm, "D_was": 98*mm, "t_was": 8.0*mm},  # 2nd choice
        "M56" : {"D":56*mm, "Pc":5.50*mm, "Pf":4.00, "t_nut":45.0*mm, "t_hex":35.0*mm, "D_hex": 93.56*mm, "d_hex": 85*mm, "d_was": 58.0*mm, "D_was":105*mm, "t_was": 9.0*mm},
        "M60" : {"D":60*mm, "Pc":5.50*mm, "Pf":4.00, "t_nut":48.0*mm, "t_hex":38.0*mm, "D_hex": 99.21*mm, "d_hex": 90*mm, "d_was": 62.0*mm, "D_was":110*mm, "t_was": 9.0*mm},  # 2nd choice
        "M64" : {"D":64*mm, "Pc":6.00*mm, "Pf":4.00, "t_nut":51.0*mm, "t_hex":40.0*mm, "D_hex":104.86*mm, "d_hex": 95*mm, "d_was": 66.0*mm, "D_was":115*mm, "t_was": 9.0*mm},
        "M72" : {"D":72*mm, "Pc":6.00*mm, "Pf":4.00, "t_nut":58.0*mm, "t_hex":46.5*mm, "D_hex":116.20*mm, "d_hex":105*mm, "d_was": 74.0*mm, "D_was":125*mm, "t_was":10.0*mm},
        "M80" : {"D":80*mm, "Pc":6.00*mm, "Pf":4.00, "t_nut":64.0*mm, "t_hex":51.6*mm, "D_hex":127.50*mm, "d_hex":115*mm, "d_was": 82.0*mm, "D_was":140*mm, "t_was":12.0*mm},
        "M90" : {"D":90*mm, "Pc":6.00*mm, "Pf":4.00, "t_nut":72.0*mm, "t_hex":57.7*mm, "D_hex":144.10*mm, "d_hex":130*mm, "d_was": 93.0*mm, "D_was":160*mm, "t_was":12.0*mm},
        "M100": {"D":100*mm,"Pc":6.00*mm, "Pf":4.00, "t_nut":80.0*mm, "t_hex":63.9*mm, "D_hex":161.02*mm, "d_hex":145*mm, "d_was":104.0*mm, "D_was":175*mm, "t_was":14.0*mm}
    },

    "materials": {

        # Carbon-steel and alloy steel bolts,
        # according to ISO 898-1

        "4.6": {
            "fy": 240*MPa,    # yield stress (stress at 0.2% non-proportional elongation)
            "fu": 400*MPa,    # ultimate tensile stress
            "E" : 210*GPa,    # elastic modulus
            "nu": 0.3         # poissons ratio
        },

        "4.8"   : {"fy": 320*MPa, "fu": 400*MPa, "E":210*GPa, "nu":0.3},
        "5.6"   : {"fy": 300*MPa, "fu": 500*MPa, "E":210*GPa, "nu":0.3},
        "5.8"   : {"fy": 400*MPa, "fu": 500*MPa, "E":210*GPa, "nu":0.3},
        "6.8"   : {"fy": 480*MPa, "fu": 600*MPa, "E":210*GPa, "nu":0.3},
        "8.8"   : {"fy": 640*MPa, "fu": 800*MPa, "E":210*GPa, "nu":0.3},
        "9.8"   : {"fy": 720*MPa, "fu": 900*MPa, "E":210*GPa, "nu":0.3},
        "10.9"  : {"fy": 900*MPa, "fu":1000*MPa, "E":210*GPa, "nu":0.3},
        "12.9"  : {"fy":1080*MPa, "fu":1200*MPa, "E":210*GPa, "nu":0.3},

        # Stainless-steel bolts,
        # according to ISO 3506-1

        # ISO 3506-1 Austenitic
        "A50" : {"fy": 210*MPa, "fu": 500*MPa, "E":210*GPa, "nu":0.3},
        "A70" : {"fy": 450*MPa, "fu": 700*MPa, "E":210*GPa, "nu":0.3},
        "A80" : {"fy": 600*MPa, "fu": 800*MPa, "E":210*GPa, "nu":0.3},
        "A100": {"fy": 800*MPa, "fu":1000*MPa, "E":210*GPa, "nu":0.3},

        # ISO 3506-1 Duplex
        "D70" : {"fy": 450*MPa, "fu": 700*MPa, "E":210*GPa, "nu":0.3},
        "D80" : {"fy": 600*MPa, "fu": 800*MPa, "E":210*GPa, "nu":0.3},
        "D100": {"fy": 800*MPa, "fu":1000*MPa, "E":210*GPa, "nu":0.3},

        # ISO 3506-1 Martensitic
        "C50" : {"fy": 250*MPa, "fu": 500*MPa, "E":210*GPa, "nu":0.3},
        "C70" : {"fy": 410*MPa, "fu": 700*MPa, "E":210*GPa, "nu":0.3},
        "C80" : {"fy": 640*MPa, "fu": 800*MPa, "E":210*GPa, "nu":0.3},
        "C110": {"fy": 820*MPa, "fu":1100*MPa, "E":210*GPa, "nu":0.3},

        # ISO 3506-1 Ferritic
        "F45" : {"fy": 250*MPa, "fu": 450*MPa, "E":210*GPa, "nu":0.3},
        "F60" : {"fy": 410*MPa, "fu": 600*MPa, "E":210*GPa, "nu":0.3}
    }
}



