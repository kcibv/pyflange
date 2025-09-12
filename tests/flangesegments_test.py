import pytest
from pyflange.flangesegments import *
from pyflange.bolts import MetricBolt, HexNut
from pyflange.gap import gap_height_distribution

from math import *
import numpy as np

# Units of measurement
deg = pi/180



class TestPolynomialLFlangeSegment:

    def fseg (self, gap_angle=30*deg, gap_shape_factor=1.0, tilt_angle=0.0):
        D = 7.5
        Nb = 120
        
        return PolynomialLFlangeSegment(

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

            gap_height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
            gap_angle = gap_angle,  # longitudinal gap length
            gap_shape_factor = gap_shape_factor,

            s_ratio = 100/72)        # ratio of bottom shell thickness over tower shell thickness


    def test_shell_force_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -123.3

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_rest/1000, 1) == -123.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -123.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -123.3


    def test_bolt_force_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2800.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2800.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2800.0


    def test_bolt_moment_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -17.2
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -21.4
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -23.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -24.3

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -17.2
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -21.4
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -23.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -24.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -17.2
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -21.4
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -23.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -24.3


    def test_shell_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 160.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 128.6
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 119.7
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 112.9

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 160.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 128.6
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 119.7
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 112.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 160.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 128.6
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 119.7
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 112.9


    def test_bolt_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2824.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2818.7
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2820.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2818.6

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2828.8
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2822.4
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2824.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2822.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2800.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2800.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2800.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2800.0


    def test_bolt_moment_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 115.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) ==  95.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 100.1
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) ==  94.7

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 134.3
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 109.6
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 115.6
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 109.2

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 22.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 22.4
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 22.4
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 22.4


    def test_shell_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1876.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1639.3
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1432.7
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1389.5

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1876.7
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1639.3
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1432.7
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1389.5

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2013.6
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2013.6
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2013.6
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2013.6


    def test_bolt_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3640.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3640.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3640.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3640.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3500.0


    def test_bolt_moment_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2986.9
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3009.7
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2992.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2999.8

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3531.9
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3554.8
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3537.3
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3544.8

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 3005.9
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 3074.6
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 3100.4
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 3123.0


    def test_shell_force_at_closed_gap (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -956.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -853.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -931.9
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -897.8

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -956.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -853.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -931.9
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -897.8

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -124.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -124.3
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -124.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -124.3


    def test_bolt_axial_force (self):

        def test (fseg, expected_Fs4):
            Z1 = fseg.shell_force_at_rest
            Fs1 = fseg.bolt_force_at_rest
            Z2 = fseg.shell_force_at_tensile_ULS
            Fs2 = fseg.bolt_force_at_tensile_ULS
            Z3 = fseg.shell_force_at_small_displacement
            Fs3 = fseg.bolt_force_at_small_displacement
            Z4 = fseg.shell_force_at_closed_gap
            Fs4 = fseg.bolt_axial_force(Z4)

            assert round(fseg.bolt_axial_force(Z1)) == round(Fs1)
            assert round(fseg.bolt_axial_force(Z2)) == round(Fs2)
            assert round(fseg.bolt_axial_force(Z3)) == round(Fs3)
            assert round(Fs4/1000, 1) == expected_Fs4

            Z = np.array([Z1, Z2, Z3])
            Fs = np.array([Fs1, Fs2, Fs3])
            assert np.all(np.abs(Fs - fseg.bolt_axial_force(Z)) < 0.1)

        test(self.fseg( 30*deg, 1.0, 0*deg), 2783.0)
        test(self.fseg( 60*deg, 1.0, 0*deg), 2792.6)
        test(self.fseg( 90*deg, 1.0, 0*deg), 2794.3)
        test(self.fseg(120*deg, 1.0, 0*deg), 2797.0)

        test(self.fseg( 30*deg, 1.2, 0*deg), 2779.6)
        test(self.fseg( 60*deg, 1.2, 0*deg), 2791.1)
        test(self.fseg( 90*deg, 1.2, 0*deg), 2793.1)
        test(self.fseg(120*deg, 1.2, 0*deg), 2796.4)

        test(self.fseg( 30*deg, 1.0, 1*deg), 2800.0)
        test(self.fseg( 60*deg, 1.0, 1*deg), 2800.0)
        test(self.fseg( 90*deg, 1.0, 1*deg), 2800.0)
        test(self.fseg(120*deg, 1.0, 1*deg), 2800.0)


    def test_bolt_bending_moment (self):

        def test (fseg, expected_Ms4):
            Z1 = fseg.shell_force_at_rest
            Ms1 = fseg.bolt_moment_at_rest
            Z2 = fseg.shell_force_at_tensile_ULS
            Ms2 = fseg.bolt_moment_at_tensile_ULS
            Z3 = fseg.shell_force_at_small_displacement
            Ms3 = fseg.bolt_moment_at_small_displacement
            Z4 = fseg.shell_force_at_closed_gap
            Ms4 = fseg.bolt_bending_moment(Z4)

            assert round(fseg.bolt_bending_moment(Z1), 1) == round(Ms1, 1)
            assert round(fseg.bolt_bending_moment(Z2), 1) == round(Ms2, 1)
            assert round(fseg.bolt_bending_moment(Z3), 1) == round(Ms3, 1)
            assert round(Ms4, 1) == expected_Ms4

            Z = np.array([Z1, Z2, Z3, Z4])
            Ms = np.array([Ms1, Ms2, Ms3, Ms4])
            assert np.all(np.abs(Ms - fseg.bolt_bending_moment(Z)) < 0.1)

        test(self.fseg( 30*deg, 1.0, 0*deg), -141.4)
        test(self.fseg( 60*deg, 1.0, 0*deg), -113.5)
        test(self.fseg( 90*deg, 1.0, 0*deg), -120.5)
        test(self.fseg(120*deg, 1.0, 0*deg), -112.4)

        test(self.fseg( 30*deg, 1.2, 0*deg), -154.7)
        test(self.fseg( 60*deg, 1.2, 0*deg), -119.3)
        test(self.fseg( 90*deg, 1.2, 0*deg), -125.0)
        test(self.fseg(120*deg, 1.2, 0*deg), -114.7)

        test(self.fseg( 30*deg, 1.0, 1*deg), -17.2)
        test(self.fseg( 60*deg, 1.0, 1*deg), -21.4)
        test(self.fseg( 90*deg, 1.0, 1*deg), -23.0)
        test(self.fseg(120*deg, 1.0, 1*deg), -24.4)


    def test_failure_mode (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fm, Zus = fseg.failure_mode(335e6, 285e6)
        assert fm == "B"



class TestPolynomialTFlangeSegment:

    def fseg (self, gap_angle=30*deg, gap_shape_factor=1.0, tilt_angle=0.0):
        D = 7.5
        Nb = 200

        M48 = MetricBolt(
            nominal_diameter = 0.048,
            thread_pitch = 0.005,
            shank_diameter_ratio = 44.752/48,
            shank_length = 0.150,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            stud = True)

        M48_hex_nut = HexNut(
            nominal_diameter = 0.048,
            thickness = 0.064,
            inscribed_diameter = 0.075,
            circumscribed_diameter = 0.0826,
            bearing_diameter = 0.092
        )

        return PolynomialTFlangeSegment(

            a = 0.0625,         # distance between inner face of the flange and center of the bolt hole
            b = 0.1110,         # distance between center of the bolt hole and center-line of the shell
            s = 0.0900,         # shell thickness
            t = 0.1200,         # flange thickness
            R = D/2,            # shell outer curvature radius
            central_angle = 2*pi/Nb,    # angle subtended by the flange segment arc

            Zg = -81400*1.3,  # load applied to the flange segment shell at rest
                          # (normally dead weight of tower + RNA, divided by the number of bolts)

            bolt = M48,
            Fv = 928000,        # applied bolt preload

            Do = 0.052,     # bolt hole diameter
            washer = None,      # no washer
            nut = M48_hex_nut,  # bolt nut

            tilt_angle = tilt_angle,

            gap_height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
            gap_angle = gap_angle,  # longitudinal gap length
            gap_shape_factor = gap_shape_factor,

            s_ratio = 1.0)        # ratio of bottom shell thickness over tower shell thickness


    def test_shell_force_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -105.8

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -105.8

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -105.8
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -105.8


    def test_bolt_force_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 928.0

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_force_at_rest/1000, 1) == 928.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 928.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 928.0


    def test_bolt_moment_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -11.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -11.4
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -12.1
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -12.0

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == -11.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == -11.4
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == -12.1
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == -12.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -14.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -14.7
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -15.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -15.2


    def test_shell_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.9
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 119.2
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 118.3
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 117.5

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.9
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 119.2
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 118.3
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 117.5

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 121.9
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 119.2
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 118.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 117.5


    def test_bolt_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 938.9        
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 939.1
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 940.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 939.7

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 944.3
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 944.7
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 946.0
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 945.6

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 943.7
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 943.9
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 944.8
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 944.5


    def test_bolt_moment_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 21.2
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 21.3
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 22.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 21.9

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 23.6
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 23.8
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 24.8
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 24.4

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 26.2
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 26.2
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 26.9
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 26.6


    def test_shell_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1633.5
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1569.7
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1491.1
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1497.6

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1633.5
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1569.7
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1491.1
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1497.6

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1274.9
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1218.8
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1143.1
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1151.9


    def test_bolt_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1276.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1276.0
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1276.0
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 1276.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 1160.0


    def test_bolt_moment_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 585.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 559.7
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 525.6
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 529.5

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 601.8
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 577.1
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 544.4
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 548.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 431.8
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 413.2
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 386.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 390.3


    def test_shell_force_at_closed_gap (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -871.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -898.8
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -963.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -946.9

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -871.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -898.8
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -963.2
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -946.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1212.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1240.1
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1304.6
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1288.2


    def test_bolt_axial_force (self):

        def test (fseg, expected_Fs4):
            Z1 = fseg.shell_force_at_rest
            Fs1 = fseg.bolt_force_at_rest
            Z2 = fseg.shell_force_at_tensile_ULS
            Fs2 = fseg.bolt_force_at_tensile_ULS
            Z3 = fseg.shell_force_at_small_displacement
            Fs3 = fseg.bolt_force_at_small_displacement
            Z4 = fseg.shell_force_at_closed_gap
            Fs4 = fseg.bolt_axial_force(Z4)

            assert round(fseg.bolt_axial_force(Z1)) == round(Fs1)
            assert round(fseg.bolt_axial_force(Z2)) == round(Fs2)
            assert round(fseg.bolt_axial_force(Z3)) == round(Fs3)
            assert round(Fs4/1000, 1) == expected_Fs4

            Z = np.array([Z1, Z2, Z3])
            Fs = np.array([Fs1, Fs2, Fs3])
            assert np.all(np.abs(Fs - fseg.bolt_axial_force(Z)) < 0.1)

        test(self.fseg( 30*deg, 1.0, 0*deg), 914.7)
        test(self.fseg( 60*deg, 1.0, 0*deg), 913.8)
        test(self.fseg( 90*deg, 1.0, 0*deg), 911.5)
        test(self.fseg(120*deg, 1.0, 0*deg), 912.2)

        test(self.fseg( 30*deg, 1.5, 0*deg), 908.0)
        test(self.fseg( 60*deg, 1.5, 0*deg), 906.8)
        test(self.fseg( 90*deg, 1.5, 0*deg), 903.2)
        test(self.fseg(120*deg, 1.5, 0*deg), 904.3)

        test(self.fseg( 30*deg, 1.0, 1*deg), 900.6)
        test(self.fseg( 60*deg, 1.0, 1*deg), 899.9)
        test(self.fseg( 90*deg, 1.0, 1*deg), 897.7)
        test(self.fseg(120*deg, 1.0, 1*deg), 898.5)


    def test_bolt_bending_moment (self):

        def test (fseg, expected_Ms4):
            Z1 = fseg.shell_force_at_rest
            Ms1 = fseg.bolt_moment_at_rest
            Z2 = fseg.shell_force_at_tensile_ULS
            Ms2 = fseg.bolt_moment_at_tensile_ULS
            Z3 = fseg.shell_force_at_small_displacement
            Ms3 = fseg.bolt_moment_at_small_displacement
            Z4 = fseg.shell_force_at_closed_gap
            Ms4 = fseg.bolt_bending_moment(Z4)

            assert round(fseg.bolt_bending_moment(Z1), 1) == round(Ms1, 1)
            assert round(fseg.bolt_bending_moment(Z2), 1) == round(Ms2, 1)
            assert round(fseg.bolt_bending_moment(Z3), 1) == round(Ms3, 1)
            assert round(Ms4, 1) == expected_Ms4

            Z = np.array([Z1, Z2, Z3, Z4])
            Ms = np.array([Ms1, Ms2, Ms3, Ms4])
            assert np.all(np.abs(Ms - fseg.bolt_bending_moment(Z)) < 0.1)

        '''
        These values do not correspond directly with the SGRE validation values.
        The deviations are due to the different calculation of the initial slope.
        The pyflange results are used as comparative values.
        '''

        test(self.fseg( 30*deg, 1.0, 0*deg), -53.3)
        test(self.fseg( 60*deg, 1.0, 0*deg), -57.0)
        test(self.fseg( 90*deg, 1.0, 0*deg), -64.9)
        test(self.fseg(120*deg, 1.0, 0*deg), -63.0)

        test(self.fseg( 30*deg, 1.5, 0*deg), -57.5)
        test(self.fseg( 60*deg, 1.5, 0*deg), -61.3)
        test(self.fseg( 90*deg, 1.5, 0*deg), -69.8)
        test(self.fseg(120*deg, 1.5, 0*deg), -67.8)

        test(self.fseg( 30*deg, 1.0, 1*deg), -96.8)
        test(self.fseg( 60*deg, 1.0, 1*deg), -101.3)
        test(self.fseg( 90*deg, 1.0, 1*deg), -110.7)
        test(self.fseg(120*deg, 1.0, 1*deg), -108.6)


    def test_failure_mode (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fm, Zus = fseg.failure_mode(335e6, 285e6)
        assert fm == "A"



def test_shell_stiffness ():
    pass
    assert shell_stiffness(1.5, 0.01, 135*deg) == 83000000



