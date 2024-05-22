
import pytest
from pyflange.flangesegments import PolynomialLFlangeSegment
from pyflange.bolts import MetricBolt
from pyflange.gap import gap_height_distribution

from math import *
import numpy as np


class TestPolynomialLFlangeSegment:

    @pytest.fixture
    def fseg (self):
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

            gap_height = gap_height_distribution(D, 0.0014, pi/6*D/2).ppf(0.95),   # maximum longitudinal gap height
            gap_angle = pi/6,  # longitudinal gap length

            s_ratio = 100/72)        # ratio of bottom shell thickness over tower shell thickness


    def test_shell_force_at_rest (self, fseg):
        assert round(fseg.shell_force_at_rest/1000, 1) == -123.3


    def test_bolt_force_at_rest (self, fseg):
        assert round(fseg.bolt_force_at_rest/1000, 1) == 2876.0


    def test_bolt_moment_at_rest (self, fseg):
        assert round(fseg.bolt_moment_at_rest, 1) == -14.7


    def test_shell_force_at_small_displacement (self, fseg):
        assert round(fseg.shell_force_at_small_displacement/1000, 1) == 251.2


    def test_bolt_force_at_small_displacement (self, fseg):
        assert round(fseg.bolt_force_at_small_displacement/1000, 1) == 2929.3


    def test_bolt_moment_at_small_displacement (self, fseg):
        assert round(fseg.bolt_moment_at_small_displacement, 1) == 195.9


    def test_shell_force_at_tensile_ULS (self, fseg):
        assert round(fseg.shell_force_at_tensile_ULS/1000, 1) == 2001.0


    def test_bolt_force_at_tensile_ULS (self, fseg):
        assert round(fseg.bolt_force_at_tensile_ULS/1000, 1) == 3595.0


    def test_bolt_moment_at_tensile_ULS (self, fseg):
        assert round(fseg.bolt_moment_at_tensile_ULS, 1) == 2477.8


    def test_shell_force_at_closed_gap (self, fseg):
        assert round(fseg.shell_force_at_closed_gap/1000, 1) == -1244.4


    def test_bolt_axial_force (self, fseg):
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
        assert round(Fs4/1000, 1) == 2814.2

        Z = np.array([Z1, Z2, Z3])
        Fs = np.array([Fs1, Fs2, Fs3])
        assert np.all(Fs - fseg.bolt_axial_force(Z) < 0.1)


    def test_bolt_bending_moment (self, fseg):
        Z1 = fseg.shell_force_at_rest
        Ms1 = fseg.bolt_moment_at_rest
        Z2 = fseg.shell_force_at_tensile_ULS
        Ms2 = fseg.bolt_moment_at_tensile_ULS
        Z3 = fseg.shell_force_at_small_displacement
        Ms3 = fseg.bolt_moment_at_small_displacement
        Z4 = fseg.shell_force_at_closed_gap
        Ms4 = fseg.bolt_bending_moment(Z4)

        assert fseg.bolt_bending_moment(Z1) == Ms1
        assert fseg.bolt_bending_moment(Z2) == Ms2
        assert fseg.bolt_bending_moment(Z3) == Ms3
        assert round(Ms4, 1) == -280.6

        Z = np.array([Z1, Z2, Z3, Z4])
        Ms = np.array([Ms1, Ms2, Ms3, Ms4])
        assert np.all(Ms == fseg.bolt_bending_moment(Z))


    def test_failure_mode (self, fseg):
        fm = fseg.failure_mode(335e6, 285e6)
        assert fm == "B"
