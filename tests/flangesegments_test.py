
import pytest
from pyflange.flangesegments import PolynomialLFlangeSegment, PolynomialTFlangeSegment
from pyflange.bolts import MetricBolt
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
            Fv = 2876000,        # applied bolt preload

            Do = 0.086,     # bolt hole diameter
            Dw = 0.140,    # washer diameter

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
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_rest/1000, 1) == 2876.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2876.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_rest/1000, 1) == 2876.0


    def test_bolt_moment_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -14.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -26.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -30.3
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -32.1

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -14.7
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -26.5
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -30.3
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -32.1

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -14.7
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -26.5
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -30.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -32.1


    def test_shell_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 251.2
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 139.1
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 114.9

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 251.2
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 139.1
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 114.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 251.2
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 139.1
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 114.9


    def test_bolt_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2929.3
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2913.3
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2910.9
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2909.9

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2939.9
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2920.8
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2917.8
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2916.7

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2893.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0


    def test_bolt_moment_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 195.9
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 146.2
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 138.5
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 135.5

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 229.1
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 169.4
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 160.2
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 156.6

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 83.9
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 29.9
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 29.9
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 29.9


    def test_shell_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 2001.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1696.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1595.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1541.2

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 2001.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1696.5
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1595.2
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1541.2

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2068.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2068.3
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2068.3
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 2068.3


    def test_bolt_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3738.8
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3738.8
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3738.8
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_tensile_ULS/1000, 1) == 3738.8

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_tensile_ULS/1000, 1) == 3595.0


    def test_bolt_moment_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2477.8
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2604.4
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2631.6
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2640.6

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2925.7
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3052.3
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3079.5
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3088.5

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2485.8
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2684.4
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2747.9
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2777.8


    def test_shell_force_at_closed_gap (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -1367.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -929.9
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -836.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -804.6

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -1367.7
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -929.9
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -836.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -804.6

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -451.4
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

        test(self.fseg( 30*deg, 1.0, 0*deg), 2813.6)
        test(self.fseg( 60*deg, 1.0, 0*deg), 2835.8)
        test(self.fseg( 90*deg, 1.0, 0*deg), 2841.6)
        test(self.fseg(120*deg, 1.0, 0*deg), 2844.0)

        test(self.fseg( 30*deg, 1.2, 0*deg), 2801.1)
        test(self.fseg( 60*deg, 1.2, 0*deg), 2827.7)
        test(self.fseg( 90*deg, 1.2, 0*deg), 2834.8)
        test(self.fseg(120*deg, 1.2, 0*deg), 2837.6)

        test(self.fseg( 30*deg, 1.0, 1*deg), 2877.9)
        test(self.fseg( 60*deg, 1.0, 1*deg), 2876.0)
        test(self.fseg( 90*deg, 1.0, 1*deg), 2876.0)
        test(self.fseg(120*deg, 1.0, 1*deg), 2876.0)


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

        test(self.fseg( 30*deg, 1.0, 0*deg), -283.2)
        test(self.fseg( 60*deg, 1.0, 0*deg), -238.5)
        test(self.fseg( 90*deg, 1.0, 0*deg), -224.9)
        test(self.fseg(120*deg, 1.0, 0*deg), -220.4)

        test(self.fseg( 30*deg, 1.2, 0*deg), -322.1)
        test(self.fseg( 60*deg, 1.2, 0*deg), -263.6)
        test(self.fseg( 90*deg, 1.2, 0*deg), -246.3)
        test(self.fseg(120*deg, 1.2, 0*deg), -240.4)

        test(self.fseg( 30*deg, 1.0, 1*deg), -28.2)
        test(self.fseg( 60*deg, 1.0, 1*deg), -26.6)
        test(self.fseg( 90*deg, 1.0, 1*deg), -30.4)
        test(self.fseg(120*deg, 1.0, 1*deg), -32.2)


    def test_failure_mode (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fm, Zus = fseg.failure_mode(335e6, 285e6)
        assert fm == "B"




class TestPolynomialTFlangeSegment:

    def fseg (self, gap_angle=30*deg, gap_shape_factor=1.0, tilt_angle=0.0):
        D = 7.5
        Nb = 200

        return PolynomialTFlangeSegment(

            a = 0.0625,         # distance between inner face of the flange and center of the bolt hole
            b = 0.1110,         # distance between center of the bolt hole and center-line of the shell
            s = 0.0900,         # shell thickness
            t = 0.1200,         # flange thickness
            R = D/2,            # shell outer curvature radius
            central_angle = 2*pi/Nb,    # angle subtended by the flange segment arc

            Zg = -81400,  # load applied to the flange segment shell at rest
                          # (normally dead weight of tower + RNA, divided by the number of bolts)

            bolt = MetricBolt(
                nominal_diameter = 0.048,
                thread_pitch = 0.005,
                shank_diameter_ratio = 44.752/48,
                shank_length = 0.150,
                yield_stress = 900e6,
                ultimate_tensile_stress = 1000e6,
                stud = True),
            Fv = 928000,        # applied bolt preload

            Do = 0.052,     # bolt hole diameter
            Dw = 0.092,    # washer diameter

            tilt_angle = tilt_angle,

            gap_height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
            gap_angle = gap_angle,  # longitudinal gap length
            gap_shape_factor = gap_shape_factor,

            s_ratio = 1.0)        # ratio of bottom shell thickness over tower shell thickness


    def test_shell_force_at_rest (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_rest/1000, 1) == -81.4

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_rest/1000, 1) == -81.4

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -81.4
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_rest/1000, 1) == -81.4


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
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == 0.0

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_rest, 1) == 0.0

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == 0.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == 0.0


    def test_shell_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 128.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 120.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 118.5
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 117.6

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 128.0
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 120.5
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 118.5
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_small_displacement/1000, 1) == 117.6

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 128.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 120.5
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 118.5
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 117.6


    def test_bolt_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 940.1
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 939.7
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 939.6
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 939.5

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 946.2
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 945.6
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 945.4
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_force_at_small_displacement/1000, 1) == 945.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 940.1
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 939.7
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 939.6
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 939.5


    def test_bolt_moment_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 14.4
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 13.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 13.3
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 13.2

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 17.9
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 16.9
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 16.6
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_small_displacement, 1) == 16.5

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 16.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 15.1
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 14.8
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 14.6


    def test_shell_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1486.2
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1539.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1554.5
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1557.9

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1486.2
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1539.0
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1554.5
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1557.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1109.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1184.3
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1205.9
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_tensile_ULS/1000, 1) == 1211.9


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
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 647.5
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 688.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 700.9
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 704.8

        assert round(self.fseg( 30*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 647.5
        assert round(self.fseg( 60*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 688.5
        assert round(self.fseg( 90*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 700.9
        assert round(self.fseg(120*deg, 1.5, 0*deg).bolt_moment_at_tensile_ULS, 1) == 704.8

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 483.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 529.8
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 543.7
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 548.3


    def test_shell_force_at_closed_gap (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -1054.8
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -920.3
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -879.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -864.4

        assert round(self.fseg( 30*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -1054.8
        assert round(self.fseg( 60*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -920.3
        assert round(self.fseg( 90*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -879.2
        assert round(self.fseg(120*deg, 1.5, 0*deg).shell_force_at_closed_gap/1000, 1) == -864.4

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1396.2
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1261.6
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1220.5
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -1205.7


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

        test(self.fseg( 30*deg, 1.0, 0*deg), 906.6)
        test(self.fseg( 60*deg, 1.0, 0*deg), 908.8)
        test(self.fseg( 90*deg, 1.0, 0*deg), 909.5)
        test(self.fseg(120*deg, 1.0, 0*deg), 909.8)

        test(self.fseg( 30*deg, 1.5, 0*deg), 895.8)
        test(self.fseg( 60*deg, 1.5, 0*deg), 899.2)
        test(self.fseg( 90*deg, 1.5, 0*deg), 900.3)
        test(self.fseg(120*deg, 1.5, 0*deg), 900.8)

        test(self.fseg( 30*deg, 1.0, 1*deg), 909.1)
        test(self.fseg( 60*deg, 1.0, 1*deg), 907.9)
        test(self.fseg( 90*deg, 1.0, 1*deg), 907.8)
        test(self.fseg(120*deg, 1.0, 1*deg), 907.8)


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

        test(self.fseg( 30*deg, 1.0, 0*deg), -7.7)
        test(self.fseg( 60*deg, 1.0, 0*deg), -6.7)
        test(self.fseg( 90*deg, 1.0, 0*deg), -6.4)
        test(self.fseg(120*deg, 1.0, 0*deg), -6.2)

        test(self.fseg( 30*deg, 1.5, 0*deg), -17.0)
        test(self.fseg( 60*deg, 1.5, 0*deg), -14.7)
        test(self.fseg( 90*deg, 1.5, 0*deg), -14.0)
        test(self.fseg(120*deg, 1.5, 0*deg), -13.8)

        test(self.fseg( 30*deg, 1.0, 1*deg), -5.7)
        test(self.fseg( 60*deg, 1.0, 1*deg), -5.7)
        test(self.fseg( 90*deg, 1.0, 1*deg), -5.6)
        test(self.fseg(120*deg, 1.0, 1*deg), -5.6)


    def test_failure_mode (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fm, Zus = fseg.failure_mode(335e6, 285e6)
        assert fm == "A"
