import pytest
from pyflange.flangesegments import *
from pyflange.bolts import MetricBolt, HexNut
from pyflange.gap import gap_height_distribution
from metrum.units import deg

from math import *
import numpy as np

import os
import sys
import json
import glob

VALIDATIONS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "validations")
sys.path.insert(0, VALIDATIONS_DIR)
from validate import validate


def validation_cases (flange_type):
    '''Paths of the validation case files in validations/references for the given flange type ("L" or "T").'''
    paths = []
    for path in sorted(glob.glob(os.path.join(VALIDATIONS_DIR, "references", "*.json"))):
        with open(path, encoding="utf-8") as f:
            if json.load(f)["segment"]["flange_type"] == flange_type:
                paths.append(path)
    return paths


def check_validation_case (path):
    failures = [c for c in validate(path) if not c["pass"]]
    assert not failures, f"{len(failures)} failed comparisons in {os.path.basename(path)}:\n" + "\n".join(
        f"    gap {c['gap_angle']} deg, {c['parameter']}: expected {c['expected']:.6g}, obtained {c['obtained']:.6g}, "
        f"mismatch {c['mismatch']:+.3%} (allowable {c['allowable_mismatch']:.2%})" for c in failures)



class TestGap:

    def test_attributes (self):
        gap = Gap(height=0.001, angle=30*deg, shape_factor=1.2)
        assert gap.height == 0.001
        assert gap.angle == 30*deg
        assert gap.shape_factor == 1.2

    def test_default_shape_factor (self):
        gap = Gap(height=0.001, angle=30*deg)
        assert gap.shape_factor == 1.0



class TestFlangeSegment:

    def test_abstract (self):
        with pytest.raises(TypeError):
            FlangeSegment()

        class IncompleteFlangeSegment (FlangeSegment):
            def bolt_axial_force (self, shell_pull):
                return shell_pull

        with pytest.raises(TypeError):
            IncompleteFlangeSegment()

    def test_concrete_subclass (self):
        class LinearFlangeSegment (FlangeSegment):
            def bolt_axial_force (self, shell_pull):
                return 2*shell_pull
            def bolt_bending_moment (self, shell_pull):
                return 3*shell_pull

        fseg = LinearFlangeSegment()
        assert fseg.bolt_axial_force(1.0) == 2.0
        assert fseg.bolt_bending_moment(1.0) == 3.0



class TestPolynomialFlangeSegment:

    class StubFlangeSegment (PolynomialFlangeSegment):
        # Reference points of the polynomial transfer functions
        shell_force_at_rest = 0.0
        bolt_force_at_rest = 100.0
        bolt_moment_at_rest = 1.0
        shell_force_at_small_displacement = 10.0
        bolt_force_at_small_displacement = 110.0
        bolt_moment_at_small_displacement = 2.0
        shell_force_at_tensile_ULS = 100.0
        bolt_force_at_tensile_ULS = 300.0
        bolt_moment_at_tensile_ULS = 30.0
        shell_force_at_closed_gap = -50.0

    def test_abstract (self):
        with pytest.raises(TypeError):
            PolynomialFlangeSegment()

    def test_bolt_axial_force (self):
        fseg = self.StubFlangeSegment()

        # The tensile branch passes through the points P1, P3 and P2
        assert fseg.bolt_axial_force(0.0) == pytest.approx(100.0)
        assert fseg.bolt_axial_force(10.0) == pytest.approx(110.0)
        assert fseg.bolt_axial_force(100.0) == pytest.approx(300.0)

        # The compressive branch is tangent to the tensile branch in P1
        # and has zero slope in P4.
        dZ = 1e-6
        tens_slope = (fseg.bolt_axial_force(dZ) - fseg.bolt_axial_force(0.0)) / dZ
        comp_slope = (fseg.bolt_axial_force(0.0) - fseg.bolt_axial_force(-dZ)) / dZ
        assert comp_slope == pytest.approx(tens_slope, rel=1e-4)
        assert (fseg.bolt_axial_force(-50.0+dZ) - fseg.bolt_axial_force(-50.0)) / dZ == pytest.approx(0.0, abs=1e-4)

        # The force is constant below P4
        Fs4 = fseg.bolt_axial_force(-50.0)
        assert Fs4 < 100.0
        assert fseg.bolt_axial_force(-60.0) == pytest.approx(Fs4)
        assert fseg.bolt_axial_force(-1e6) == pytest.approx(Fs4)

        # Array input
        Z = np.array([-1e6, -50.0, 0.0, 10.0, 100.0])
        np.testing.assert_allclose(fseg.bolt_axial_force(Z), [Fs4, Fs4, 100.0, 110.0, 300.0])

    def test_bolt_bending_moment (self):
        fseg = self.StubFlangeSegment()

        # The tensile branch passes through the points P1, P3 and P2
        assert fseg.bolt_bending_moment(0.0) == pytest.approx(1.0)
        assert fseg.bolt_bending_moment(10.0) == pytest.approx(2.0)
        assert fseg.bolt_bending_moment(100.0) == pytest.approx(30.0)

        # The compressive branch is tangent to the tensile branch in P1
        # and has zero slope in P4.
        dZ = 1e-6
        tens_slope = (fseg.bolt_bending_moment(dZ) - fseg.bolt_bending_moment(0.0)) / dZ
        comp_slope = (fseg.bolt_bending_moment(0.0) - fseg.bolt_bending_moment(-dZ)) / dZ
        assert comp_slope == pytest.approx(tens_slope, rel=1e-4)
        assert (fseg.bolt_bending_moment(-50.0+dZ) - fseg.bolt_bending_moment(-50.0)) / dZ == pytest.approx(0.0, abs=1e-4)

        # The moment is constant below P4
        Ms4 = fseg.bolt_bending_moment(-50.0)
        assert fseg.bolt_bending_moment(-1e6) == pytest.approx(Ms4)

        # Array input
        Z = np.array([-1e6, -50.0, 0.0, 10.0, 100.0])
        np.testing.assert_allclose(fseg.bolt_bending_moment(Z), [Ms4, Ms4, 1.0, 2.0, 30.0])



class TestPolynomialLFlangeSegment:

    def fseg (self, gap_angle=30*deg, gap_shape_factor=1.0, tilt_angle=0.0, PPL=0.05, DFT=300e-6, settlfac=0.50):
        from metrum.units import m, mm, um, kN, MPa

        D = 7500*mm
        Nb = 120
        
        return PolynomialLFlangeSegment(

            a = 232.5*mm,         # distance between inner face of the flange and center of the bolt hole
            b = 166.5*mm,         # distance between center of the bolt hole and center-line of the shell
            s =  72.0*mm,         # shell thickness
            t = 200.0*mm,         # flange thickness
            R = D/2,            # shell outer curvature radius
            central_angle = 2*pi/Nb,    # angle subtended by the flange segment arc

            Zg = -14795*kN / Nb,  # load applied to the flange segment shell at rest
                                  # (normally dead weight of tower + RNA, divided by the number of bolts)

            bolt = MetricBolt(
                nominal_diameter = 80*mm,
                thread_pitch = 6*mm,
                shank_diameter_ratio = 76.1/80,
                shank_length = 270*mm,
                yield_stress = 900*MPa,
                ultimate_tensile_stress = 1000*MPa,
                stud = True),
            Fv = 3007.85395*kN,        # applied bolt preload

            Do = 86*mm,     # bolt hole diameter
            washer = None,    # no washer diameter
            nut = HexNut(
                nominal_diameter = 80*mm,
                thickness = 64*mm,
                inscribed_diameter = 115*mm,
                circumscribed_diameter = 127.5*mm,
                bearing_diameter = 140*mm),

            tilt_angle = tilt_angle,
            
            gap = Gap(
                      height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
                      angle = gap_angle,  # longitudinal gap length
                      shape_factor = gap_shape_factor),

            s_ratio = 100/72,        # ratio of bottom shell thickness over tower shell thickness

            PPL = PPL,                 # Percentage of plastic loss
            DFT = DFT,               # Dry film thickness
            settlement_factor = settlfac)    # Bolt settlement factor

    def test_bolt_preload (self):
        from metrum.units import um, kN
        fseg = self.fseg( 30*deg, 1.0, 0*deg, PPL=0.05, DFT=300*um, settlfac=0.50)
        assert round(fseg.bolt_preload) == 2800*kN

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
        from scipy.optimize import fsolve

        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fy_sh, fy_fl, gamma_0 = 335e6, 285e6, 1.1
        fm, (Zu_A, Zu_B, Zu_D, Zu_E) = fseg.failure_mode(fy_sh, fy_fl, gamma_0)
        assert fm == "B"

        # Independent check of failure mode D (IEC 61400-6, eqs. G.4, G.7, G.9 and G.11),
        # with the segment width c at the middle surface of the shell and c' = c - dB.
        fd_sh, fd_fl = fy_sh / gamma_0, fy_fl / gamma_0
        F_tR = fseg.bolt.ultimate_tensile_capacity()
        Dw = fseg.nut.bearing_diameter
        c = fseg.central_angle * (fseg.R - fseg.s/2)
        M_pl2_hole = fd_fl * (c - fseg.Do) * fseg.t**2 / 4
        DM_pl2 = F_tR/2 * (Dw + fseg.Do)/4
        M_pl_sh = lambda N: fd_sh * c * fseg.s**2 / 4 * (1 - (N / (fd_sh * c * fseg.s))**2)
        M_pl_fl = lambda V: fd_fl * c * fseg.t**2 / 4 * (1 - (V / (fd_fl * c * fseg.t / 3**0.5))**2)**0.5
        b_fl = fseg.b - fseg.s/2 - 0.8*fseg.r
        Zu_D_sh = fsolve(lambda Z: (M_pl2_hole + DM_pl2 + M_pl_sh(Z)) / fseg.b - Z, F_tR)[0]
        Zu_D_fl = fsolve(lambda Z: (M_pl2_hole + DM_pl2 + M_pl_fl(Z)) / b_fl - Z, F_tR)[0]
        assert Zu_A == pytest.approx(F_tR)
        assert Zu_D == pytest.approx(min(Zu_D_sh, Zu_D_fl))


    def test_prying_lever_ratio (self):
        from dataclasses import replace
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)

        # Within the validity range of IEC 61400-6, eqs. G.16 and G.17, a' < a
        a, b, t = fseg.a, fseg.b, fseg.t
        alpha = t / (a + b)
        beta = (a/b - 1.25)**0.32 + 0.45
        a_red = (1 - (1 - alpha**beta)**5) * a
        assert fseg._a_red == pytest.approx(a_red)
        assert fseg._a_red < a

        # Rigid flange (a/b <= 1.25): a' = a
        assert replace(fseg, a=1.2*b)._a_red == pytest.approx(1.2*b)

        # a/b > 2.25: out of the validity range, a' = a
        assert replace(fseg, a=2.5*b)._a_red == pytest.approx(2.5*b)

        # Flange aspect ratio below -0.12*a/b + 0.55: out of the validity range, a' = a
        thin = replace(fseg, t=0.3*(a + b))
        assert thin.t / (a + b) < 0.55 - 0.12*a/b
        assert thin._a_red == pytest.approx(a)

        # Flange aspect ratio above 1: out of the validity range, a' = a
        assert replace(fseg, t=1.1*(a + b))._a_red == pytest.approx(a)


    def test_validate (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fseg.validate(335e6, 285e6)     # failure mode B: no error

        fseg.failure_mode = lambda *args, **kwargs: ("A", [])
        with pytest.raises(ValueError):
            fseg.validate(335e6, 285e6)


    @pytest.mark.parametrize("path", validation_cases("L"), ids=os.path.basename)
    def test_validation_cases (self, path):
        check_validation_case(path)



class TestPolynomialLFlangeSegmentED2:

    def fseg (self, gap_angle=30*deg, **kwargs):
        from dataclasses import fields
        fseg = TestPolynomialLFlangeSegment().fseg(gap_angle, **kwargs)
        return PolynomialLFlangeSegmentED2(**{f.name: getattr(fseg, f.name) for f in fields(fseg)})

    def test_inheritance (self):
        fseg = self.fseg()
        assert isinstance(fseg, PolynomialLFlangeSegment)

    def test_a_star (self):
        fseg = self.fseg()
        a_red = fseg._a_red
        assert fseg._a_star == pytest.approx(max(0.3, min((fseg.t / (a_red + fseg.b))**2, 1)) * a_red)

    def test_gap_stiffness (self):
        for gap_angle in (20*deg, 30*deg, 60*deg, 120*deg):
            fseg = self.fseg(gap_angle)
            f_tot = min(1.4 + 0.6 * gap_angle/(pi/2), 2.0) * min(gap_angle/(pi/6), 1)**2
            assert fseg._gap_stiffness == pytest.approx(f_tot * fseg._k_shell + fseg._k_flange)

    def test_parallel_gap_neutralization_shell_force (self):
        for gap_angle in (30*deg, 60*deg, 120*deg):
            fseg = self.fseg(gap_angle)
            s_avg = fseg.s * (1 + fseg.s_ratio) / 2
            c = fseg.central_angle * (fseg.R - s_avg/2)
            base_force = -0.5 * fseg._gap_stiffness * fseg.gap.height * c

            # The ED2 correction factor is between 1.0 and 1.5
            f_kgap = fseg._parallel_gap_neutralization_shell_force / base_force
            assert 1.0 <= f_kgap <= 1.5

    def test_bolt_axial_force (self):
        fseg = self.fseg()
        assert fseg.bolt_axial_force(fseg.shell_force_at_rest) == pytest.approx(fseg.bolt_force_at_rest)
        assert fseg.bolt_axial_force(fseg.shell_force_at_small_displacement) == pytest.approx(fseg.bolt_force_at_small_displacement)
        assert fseg.bolt_axial_force(fseg.shell_force_at_tensile_ULS) == pytest.approx(fseg.bolt_force_at_tensile_ULS)

    def test_bolt_bending_moment (self):
        fseg = self.fseg()
        assert fseg.bolt_bending_moment(fseg.shell_force_at_rest) == pytest.approx(fseg.bolt_moment_at_rest)
        assert fseg.bolt_bending_moment(fseg.shell_force_at_small_displacement) == pytest.approx(fseg.bolt_moment_at_small_displacement)
        assert fseg.bolt_bending_moment(fseg.shell_force_at_tensile_ULS) == pytest.approx(fseg.bolt_moment_at_tensile_ULS)

    def test_shell_force_at_closed_gap (self):
        # The ED2 gap closing force differs from the ED1 one
        fseg_ed1 = TestPolynomialLFlangeSegment().fseg(30*deg)
        fseg_ed2 = self.fseg(30*deg)
        assert fseg_ed2.shell_force_at_closed_gap != pytest.approx(fseg_ed1.shell_force_at_closed_gap)



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

            gap = Gap( height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
                       angle = gap_angle,  # longitudinal gap length
                       shape_factor = gap_shape_factor),

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
        The deviations are due to the different calculation of the initial slope
        and to the tensile moment curve, which PyFlange evaluates as the cubic
        polynomial of IEC 61400-6, eq. G.61, while SGRE fits a quadratic polynomial.
        The pyflange results are used as comparative values.
        '''

        test(self.fseg( 30*deg, 1.0, 0*deg), -51.6)
        test(self.fseg( 60*deg, 1.0, 0*deg), -55.1)
        test(self.fseg( 90*deg, 1.0, 0*deg), -62.9)
        test(self.fseg(120*deg, 1.0, 0*deg), -61.1)

        test(self.fseg( 30*deg, 1.5, 0*deg), -54.9)
        test(self.fseg( 60*deg, 1.5, 0*deg), -58.6)
        test(self.fseg( 90*deg, 1.5, 0*deg), -66.8)
        test(self.fseg(120*deg, 1.5, 0*deg), -64.8)

        test(self.fseg( 30*deg, 1.0, 1*deg), -94.1)
        test(self.fseg( 60*deg, 1.0, 1*deg), -98.4)
        test(self.fseg( 90*deg, 1.0, 1*deg), -107.5)
        test(self.fseg(120*deg, 1.0, 1*deg), -105.4)


    def test_tensile_moment_polynomial (self):
        # The tensile moment curve is the cubic polynomial of IEC 61400-6, eq. G.61
        for tilt_angle in (0*deg, 1*deg):
            fseg = self.fseg(60*deg, 1.0, tilt_angle)
            poly = fseg._tensile_moment_polynomial
            assert poly.degree() == 3
            Z1 = fseg.shell_force_at_rest
            Z2 = fseg.shell_force_at_tensile_ULS
            for Z in np.linspace(Z1, Z2, 7):
                expected = fseg._bolt_moment(Z, fseg.bolt_axial_force(Z))
                assert fseg.bolt_bending_moment(Z) == pytest.approx(expected)


    def test_failure_mode (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fy_fl, gamma_0 = 285e6, 1.1
        fm, (Zu_A, Zu_B, Zu_D, Zu_E) = fseg.failure_mode(335e6, fy_fl, gamma_0)
        assert fm == "A"

        # Independent check of IEC 61400-6, eqs. G.24 to G.29: each ultimate force
        # is a fixed point of the corresponding equation, with the shear force
        # in the flange equal to half of the segment force.
        fd = fy_fl / gamma_0
        F_tR = fseg.bolt.ultimate_tensile_capacity()
        Dw = fseg.nut.bearing_diameter
        c = fseg.central_angle * (fseg.R - fseg.s/2)
        b_BD = fseg.b - fseg.s/2 - 4/5*fseg.r
        b_E = fseg.b - (Dw + fseg.Do)/4 - fseg.s/2 - 4/5*fseg.r
        M_pl2_hole = fd * (c - fseg.Do) * fseg.t**2 / 4
        M_pl2 = fd * c * fseg.t**2 / 4
        DM_pl2 = F_tR/2 * (Dw + fseg.Do)/4
        V_pl = fd * fseg.t * c / 3**0.5
        M_pl3 = lambda Fu: M_pl2 * (1 - (Fu/2 / V_pl)**2)**0.5 if Fu/2 < V_pl else 0
        a = fseg._a_red

        assert Zu_A == pytest.approx(2 * F_tR)
        assert Zu_B == pytest.approx(2 * (F_tR*a + M_pl3(Zu_B)) / (a + b_BD))
        assert Zu_D == pytest.approx(2 * (M_pl2_hole + DM_pl2 + M_pl3(Zu_D)) / b_BD)
        assert Zu_E == pytest.approx(2 * (M_pl2 + M_pl3(Zu_E)) / b_E)


    def test_validate (self):
        fseg = self.fseg(30*deg, 1.0, 0.0*deg)
        fseg.validate(335e6, 285e6)     # failure mode A: no error

        fseg.failure_mode = lambda *args, **kwargs: ("B", [])
        with pytest.raises(ValueError):
            fseg.validate(335e6, 285e6)


    def test_polynomial_initial_slope (self):
        from dataclasses import replace

        def expected_slope (fseg):
            Ks = fseg._bolt_axial_stiffness
            Kp = fseg._flange_axial_stiffness
            p = Ks / (Ks + Kp)
            p_max = (fseg.bolt_force_at_tensile_ULS - fseg.bolt_force_at_rest) / (fseg.shell_force_at_tensile_ULS - fseg.shell_force_at_rest)
            chi = 0.5 * p * (-fseg._total_gap_neutralization_shell_force) / (2 * fseg.bolt_preload * 0.5)
            return min(chi, p_max) * fseg.gap.shape_factor, chi, p_max

        # Ordinary gap: IEC 61400-6, eq. G.54, below its upper bound
        fseg = self.fseg(60*deg, 1.5, 0*deg)
        slope, chi, p_max = expected_slope(fseg)
        assert chi < p_max
        assert fseg._polynomial_initial_slope == pytest.approx(slope)

        # Very large (unrealistic) gap: the slope is limited to (Fs2 - Fs1) / (Z2 - Z1)
        fseg = replace(fseg, gap=Gap(height=0.050, angle=60*deg, shape_factor=1.5))
        slope, chi, p_max = expected_slope(fseg)
        assert chi > p_max
        assert fseg._polynomial_initial_slope == pytest.approx(p_max * 1.5)


    def test_bolt_preload (self):
        from dataclasses import replace

        # No losses by default
        fseg = self.fseg(30*deg, 1.0, 0*deg)
        assert fseg.bolt_preload == fseg.Fv
        assert fseg.bolt_force_at_rest == fseg.Fv

        # Same preload losses as the L-flange
        fseg = replace(fseg, PPL=0.05, DFT=300e-6, settlement_factor=0.5)
        Ks = fseg._bolt_axial_stiffness
        Kp = fseg._flange_axial_stiffness
        f_z = max(8 * 300**(-0.6) * 300e-6, 50e-6)
        expected = fseg.Fv - 0.5 * f_z / (1/Ks + 1/Kp) - 0.05 * fseg.Fv
        assert fseg.bolt_preload == pytest.approx(expected)
        assert fseg.bolt_force_at_rest == pytest.approx(expected)
        assert fseg.bolt_axial_force(fseg.shell_force_at_rest) == pytest.approx(expected)


    def test_tilt_angle_required (self):
        from dataclasses import fields
        fseg = self.fseg(30*deg, 1.0, 0*deg)
        params = {f.name: getattr(fseg, f.name) for f in fields(fseg) if f.name != "tilt_angle"}
        with pytest.raises(TypeError):
            PolynomialTFlangeSegment(**params)


    def test_poisson_ratio (self):
        # The tilt neutralization force depends on the Poisson ratio derived from E and G
        from dataclasses import replace
        fseg = self.fseg(30*deg, 1.0, 1*deg)
        fseg_nu = replace(fseg, G=fseg.E / (2*(1 + 0.2)))
        assert fseg_nu._tilt_neutralization_shell_force != pytest.approx(fseg._tilt_neutralization_shell_force)


    def test_k_shell (self):
        from dataclasses import replace
        fseg = self.fseg(30*deg, 1.0, 0*deg)
        assert replace(fseg, k_shell=1.5e9)._k_shell == 1.5e9
        assert replace(fseg, k_shell=lambda fs: 2.5e9)._k_shell == 2.5e9

        # The default is the simplified formula of IEC 61400-6, eqs. G.43 and G.44
        L_gap = fseg.R * fseg.gap.angle
        k_fac = max(1.8, 1.3 + (8.0e-4 - 1.6e-7 * fseg.R*1000) * L_gap*1000)
        assert fseg._k_shell == pytest.approx(fseg.E * fseg.s / (k_fac * L_gap))


    @pytest.mark.parametrize("path", validation_cases("T"), ids=os.path.basename)
    def test_validation_cases (self, path):
        check_validation_case(path)



class TestPolynomialTFlangeSegmentED2:

    def fseg (self, gap_angle=30*deg, *args):
        from dataclasses import fields
        fseg = TestPolynomialTFlangeSegment().fseg(gap_angle, *args)
        return PolynomialTFlangeSegmentED2(**{f.name: getattr(fseg, f.name) for f in fields(fseg)})

    def test_inheritance (self):
        fseg = self.fseg()
        assert isinstance(fseg, PolynomialTFlangeSegment)

    def test_gap_stiffness (self):
        for gap_angle in (20*deg, 30*deg, 60*deg, 120*deg, 180*deg):
            fseg = self.fseg(gap_angle)
            f_tot = min(1.4 + 0.6 * gap_angle/(pi/2), 2.0)
            assert fseg._gap_stiffness == pytest.approx(f_tot * fseg._k_shell + fseg._k_flange)

    def test_early_prying_uses_initial_stiffness (self):
        # The edge contact preload (IEC 61400-6, eq. G.48) uses the initial linear
        # stiffness, which is the same for AMD1 and ED2.
        fseg_amd1 = TestPolynomialTFlangeSegment().fseg(60*deg)
        fseg_ed2 = self.fseg(60*deg)
        assert fseg_ed2._early_prying_neutralization_shell_force == pytest.approx(fseg_amd1._early_prying_neutralization_shell_force)

    def test_bolt_axial_force (self):
        fseg = self.fseg()
        assert fseg.bolt_axial_force(fseg.shell_force_at_rest) == pytest.approx(fseg.bolt_force_at_rest)
        assert fseg.bolt_axial_force(fseg.shell_force_at_small_displacement) == pytest.approx(fseg.bolt_force_at_small_displacement)
        assert fseg.bolt_axial_force(fseg.shell_force_at_tensile_ULS) == pytest.approx(fseg.bolt_force_at_tensile_ULS)

    def test_bolt_bending_moment (self):
        fseg = self.fseg()
        assert fseg.bolt_bending_moment(fseg.shell_force_at_rest) == pytest.approx(fseg.bolt_moment_at_rest)
        assert fseg.bolt_bending_moment(fseg.shell_force_at_small_displacement) == pytest.approx(fseg.bolt_moment_at_small_displacement)
        assert fseg.bolt_bending_moment(fseg.shell_force_at_tensile_ULS) == pytest.approx(fseg.bolt_moment_at_tensile_ULS)

    def test_shell_force_at_closed_gap (self):
        # The ED2 gap closing force differs from the AMD1 one
        fseg_amd1 = TestPolynomialTFlangeSegment().fseg(30*deg)
        fseg_ed2 = self.fseg(30*deg)
        assert fseg_ed2.shell_force_at_closed_gap != pytest.approx(fseg_amd1.shell_force_at_closed_gap)



def test_bolt_markov_matrix ():
    from pyflange.fatigue import MarkovMatrix

    # Create a test flange segment
    D = 7.5
    Nb = 120
    gap_angle = 30*deg

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

        tilt_angle = 0.0,

        gap = Gap(height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),   # maximum longitudinal gap height
                  angle = gap_angle,  # longitudinal gap length
                  shape_factor = 1.0),

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

    # Bending factor and stress multiplication factor
    mkv_bolt = bolt_markov_matrix(fseg, mkv_flng, bending_factor=0.5, SMF=1.2)
    S_min = fseg.bolt_axial_force(Z_min)/bolt_A + 0.5*fseg.bolt_bending_moment(Z_min)/bolt_W
    S_max = fseg.bolt_axial_force(Z_max)/bolt_A + 0.5*fseg.bolt_bending_moment(Z_max)/bolt_W
    np.testing.assert_allclose(mkv_bolt.mean, (S_min + S_max)/2)
    np.testing.assert_allclose(mkv_bolt.range, 1.2*np.abs(S_max - S_min))

    # Mean, range and macro-geometric factors
    mkv_bolt = bolt_markov_matrix(fseg, mkv_flng, macro_geometric_factor=1.1, mean_factor=0.9, range_factor=1.3)
    Z_mean   = mkv_flng.mean / flange_W * shell_A * 0.9 * 1.1
    Z_ranges = mkv_flng.range / flange_W * shell_A * 1.3 * 1.1
    Z_min = Z_mean - Z_ranges/2 + fseg.Zg * 1.1
    Z_max = Z_mean + Z_ranges/2 + fseg.Zg * 1.1
    S_min = fseg.bolt_axial_force(Z_min)/bolt_A
    S_max = fseg.bolt_axial_force(Z_max)/bolt_A
    np.testing.assert_allclose(mkv_bolt.mean, (S_min + S_max)/2)
    np.testing.assert_allclose(mkv_bolt.range, np.abs(S_max - S_min))



def test_shell_stiffness ():
    assert shell_stiffness(1.5, 0.01, 135*deg) == 83000000

    # No extrapolation outside the dataset range
    assert isnan(shell_stiffness(1000.0, 0.01, 135*deg))



def test_load_shell_stiffness_interpolator ():
    from pyflange.flangesegments import _load_shell_stiffness_interpolator
    from scipy.interpolate import LinearNDInterpolator

    interp = _load_shell_stiffness_interpolator()
    assert isinstance(interp, LinearNDInterpolator)

    # The interpolator is cached
    assert _load_shell_stiffness_interpolator() is interp

    assert interp(1.5, 0.01, 135*deg) == shell_stiffness(1.5, 0.01, 135*deg)



