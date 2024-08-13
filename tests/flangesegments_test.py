
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
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -14.4
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -25.4
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -29.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_rest, 1) == -30.7

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -14.4
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -25.4
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -29.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_rest, 1) == -30.7

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -14.4
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -25.4
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -29.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_rest, 1) == -30.7


    def test_shell_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 245.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 138.8
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_small_displacement/1000, 1) == 114.9

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 245.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 138.8
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_small_displacement/1000, 1) == 114.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 245.0
        assert round(self.fseg( 60*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 138.8
        assert round(self.fseg( 90*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 121.7
        assert round(self.fseg(120*deg, 1.0, 1*deg).shell_force_at_small_displacement/1000, 1) == 114.9


    def test_bolt_force_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2920.5
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2896.8
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2893.2
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2892.1

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2929.4
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2901.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2896.7
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_force_at_small_displacement/1000, 1) == 2895.3

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2887.1
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_force_at_small_displacement/1000, 1) == 2876.0


    def test_bolt_moment_at_small_displacement (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) == 167.3
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) ==  93.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) ==  82.3
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_small_displacement, 1) ==  78.7

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 195.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) == 106.4
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) ==  93.1
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_small_displacement, 1) ==  88.7

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 63.3
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 28.7
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 28.7
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_small_displacement, 1) == 28.7


    def test_shell_force_at_tensile_ULS (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 2006.0
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1698.5
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1596.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1541.0

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 2006.0
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1698.5
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1596.0
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_tensile_ULS/1000, 1) == 1541.0

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
        assert round(self.fseg( 30*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2473.9
        assert round(self.fseg( 60*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2589.8
        assert round(self.fseg( 90*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2615.0
        assert round(self.fseg(120*deg, 1.0, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2623.3

        assert round(self.fseg( 30*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 2921.8
        assert round(self.fseg( 60*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3037.7
        assert round(self.fseg( 90*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3062.9
        assert round(self.fseg(120*deg, 1.2, 0*deg).bolt_moment_at_tensile_ULS, 1) == 3071.2

        assert round(self.fseg( 30*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2481.2
        assert round(self.fseg( 60*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2666.1
        assert round(self.fseg( 90*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2726.1
        assert round(self.fseg(120*deg, 1.0, 1*deg).bolt_moment_at_tensile_ULS, 1) == 2754.6


    def test_shell_force_at_closed_gap (self):
        assert round(self.fseg( 30*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -1344.7
        assert round(self.fseg( 60*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -926.0
        assert round(self.fseg( 90*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -834.8
        assert round(self.fseg(120*deg, 1.0, 0*deg).shell_force_at_closed_gap/1000, 1) == -804.9

        assert round(self.fseg( 30*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -1344.7
        assert round(self.fseg( 60*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -926.0
        assert round(self.fseg( 90*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -834.8
        assert round(self.fseg(120*deg, 1.2, 0*deg).shell_force_at_closed_gap/1000, 1) == -804.9

        assert round(self.fseg( 30*deg, 1.0, 1*deg).shell_force_at_closed_gap/1000, 1) == -428.4
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

        test(self.fseg( 30*deg, 1.0, 0*deg), 2829.9)
        test(self.fseg( 60*deg, 1.0, 0*deg), 2865.4)
        test(self.fseg( 90*deg, 1.0, 0*deg), 2871.5)
        test(self.fseg(120*deg, 1.0, 0*deg), 2873.8)

        test(self.fseg( 30*deg, 1.2, 0*deg), 2820.7)
        test(self.fseg( 60*deg, 1.2, 0*deg), 2863.3)
        test(self.fseg( 90*deg, 1.2, 0*deg), 2870.6)
        test(self.fseg(120*deg, 1.2, 0*deg), 2873.3)

        test(self.fseg( 30*deg, 1.0, 1*deg), 2880.6)
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

        test(self.fseg( 30*deg, 1.0, 0*deg), -229.4)
        test(self.fseg( 60*deg, 1.0, 0*deg), -141.3)
        test(self.fseg( 90*deg, 1.0, 0*deg), -126.7)
        test(self.fseg(120*deg, 1.0, 0*deg), -122.5)

        test(self.fseg( 30*deg, 1.2, 0*deg), -258.2)
        test(self.fseg( 60*deg, 1.2, 0*deg), -147.9)
        test(self.fseg( 90*deg, 1.2, 0*deg), -129.5)
        test(self.fseg(120*deg, 1.2, 0*deg), -123.9)

        test(self.fseg( 30*deg, 1.0, 1*deg), -18.0)
        test(self.fseg( 60*deg, 1.0, 1*deg), -25.5)
        test(self.fseg( 90*deg, 1.0, 1*deg), -29.1)
        test(self.fseg(120*deg, 1.0, 1*deg), -30.8)


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
