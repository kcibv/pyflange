import pytest
from pyflange.fatigue import *
from pyflange.flangesegments import *
from pyflange.bolts import *
from pyflange.gap import *
from math import pi
import numpy as np
import os

# Units of measurement
m = 1
mm = 0.001*m

kg = 1
t = 1000*kg

s = 1

N = kg*m/s**2
kN = 1000*N

Pa = 1
MPa = 1e6*Pa
GPa = 1e9*Pa

rad = 1
deg = (pi/180)*rad



def test_bolt_markov_matrix ():

    # Create a test flange segment
    D = 7.5
    Nb = 120
    gap_angle=30*deg
    gap_shape_factor=1.0
    tilt_angle=0.0
    
    fseg = PolynomialLFlangeSegment(

        a = 0.2325,         # distance between inner face of the flange and center of the bolt hole
        b = 0.1665,         # distance between center of the bolt hole and center-line of the shell
        s = 0.0720,         # shell thickness
        t = 0.2000,         # flange thickness
        R = D/2,            # shell outer curvature radius
        central_angle = 2*pi/Nb,    # angle subtended by the flange segment arc

        Zg = -14795000/Nb,  # load applied to the flange segment shell at rest
                            # (normally dead weight of tower + RNA, divided by the number of bolts)

        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            shank_diameter_ratio = 76.1/80,
            shank_length = 0.270,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            stud = True),
        Fv = 2800000,        # applied bolt preload

        Do = 0.086,     # bolt hole diameter
        washer = None,    # no washer diameter
        nut = HexNut(
            nominal_diameter = 0.080,
            thickness = 0.064,
            inscribed_diameter = 0.115,
            circumscribed_diameter = 0.1275,
            bearing_diameter = 0.140),

        tilt_angle = tilt_angle,

        gap = Gap(height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
                  angle = gap_angle,  # longitudinal gap length
                  shape_factor = gap_shape_factor),

        s_ratio = 100/72)        # ratio of bottom shell thickness over tower shell thickness

    #Flange Geometry
    Rm = fseg.R - fseg.s/2
    flange_W = pi/4 * (fseg.R**4 - (fseg.R-fseg.s)**4) / Rm
    shell_A = fseg.s * fseg.central_angle * Rm

    # Bolt Geometry
    bolt_A = fseg.bolt.thread_cross_section.area
    bolt_W = fseg.bolt.thread_cross_section.elastic_section_modulus

    # Create a test flange makov matrix
    mkv_flng = MarkovMatrix(
            range = np.array([100, 200, 300]),
            mean = np.array([400, 500 ,600]),
            cycles = np.array([0.1, 0.2, 0.3])*1e6,
            duration=1)

    # Create the bolt markov matrix
    mkv_bolt = bolt_markov_matrix(fseg, mkv_flng)

    # Verify that the bolt cycles are the same as the flange cycles
    np.testing.assert_array_equal(mkv_flng.cycles, mkv_bolt.cycles)

    # Verify that the bolt load duration is the same as the flange load duration
    assert mkv_bolt.duration == mkv_flng.duration

    # Markov matrix of Z
    Z_mean   = mkv_flng.mean / flange_W * shell_A
    Z_ranges = mkv_flng.range / flange_W * shell_A
    assert round(Z_ranges[0], 3) == 0.449
    assert round(Z_ranges[1], 3) == 0.897
    assert round(Z_ranges[2], 3) == 1.346
    assert round(Z_mean[0], 3) == 1.795
    assert round(Z_mean[1], 3) == 2.244
    assert round(Z_mean[2], 3) == 2.692
    Z_min = Z_mean - Z_ranges/2 + fseg.Zg
    Z_max = Z_mean + Z_ranges/2 + fseg.Zg

    # Bolt markov matrix
    S_min = fseg.bolt_axial_force(Z_min)/bolt_A
    S_max = fseg.bolt_axial_force(Z_max)/bolt_A
    np.testing.assert_array_equal(mkv_bolt.mean, (S_min + S_max)/2)
    np.testing.assert_array_equal(mkv_bolt.range, S_max - S_min)




