import pytest
from pyflange.fatigue import *
from pyflange.flangesegments import *
from pyflange.bolts import *
from pyflange.gap import *
from math import pi, log
import numpy as np

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



def lflange_segment (gap_angle=30*deg):
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

        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            shank_diameter_ratio = 76.1/80,
            shank_length = 0.270,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            stud = True),
        Fv = 2800000,       # applied bolt preload

        Do = 0.086,         # bolt hole diameter
        washer = None,      # no washer
        nut = HexNut(
            nominal_diameter = 0.080,
            thickness = 0.064,
            inscribed_diameter = 0.115,
            circumscribed_diameter = 0.1275,
            bearing_diameter = 0.140),

        gap = Gap(height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95),
                  angle = gap_angle),

        s_ratio = 100/72)   # ratio of bottom shell thickness over tower shell thickness


def flange_markov_matrix ():
    # Bending moments on the flange, in N*m
    return MarkovMatrix(
        range  = np.array([10, 30, 60, 90]) * 1e6,
        mean   = np.array([20, 20, 30, 30]) * 1e6,
        cycles = np.array([1e7, 1e6, 1e5, 1e3]),
        duration = 20)



class TestMarkovMatrix:

    def test_attributes (self):
        mkvm = MarkovMatrix(
            range = np.array([1.0, 2.0]),
            mean = np.array([3.0, 4.0]),
            cycles = np.array([10, 20]),
            duration = 5)

        np.testing.assert_array_equal(mkvm.range, [1.0, 2.0])
        np.testing.assert_array_equal(mkvm.mean, [3.0, 4.0])
        np.testing.assert_array_equal(mkvm.cycles, [10, 20])
        assert mkvm.duration == 5

    def test_default_duration (self):
        mkvm = MarkovMatrix(range=np.array([1.0]), mean=np.array([0.0]), cycles=np.array([1]))
        assert mkvm.duration == 1.0



class TestFatigueCurve:

    class InverseSquareCurve (FatigueCurve):
        def N (self, DS):
            return 1e12 / DS**2

    def test_N (self):
        assert FatigueCurve().N(100) is None

    def test_DS (self):
        assert FatigueCurve().DS(1e6) is None

    def test_damage (self):
        curve = self.InverseSquareCurve()
        assert curve.damage(1e4, 1e3) == 1e4 / 1e6
        np.testing.assert_allclose(curve.damage(np.array([1e4, 1e2]), np.array([1e3, 1e4])), [1e-2, 1e-2])

    def test_cumulated_damage (self):
        curve = self.InverseSquareCurve()
        mkvm = MarkovMatrix(
            range = np.array([1e3, 1e4, 1e4]),
            mean = np.zeros(3),
            cycles = np.array([1e4, 1e2, np.nan]))      # NaN damages are ignored
        D = curve.cumulated_damage(mkvm)
        assert isinstance(D, float)
        assert round(D, 12) == 0.02



class TestSingleSlopeFatigueCurve:

    def test_a (self):
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6)
        assert curve.a == 2e12

    def test_N (self):
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6)
        assert curve.N(100) == 2e6
        assert curve.N(50) == 1.6e7
        np.testing.assert_allclose(curve.N(np.array([50., 100., 200.])), [1.6e7, 2e6, 2.5e5])

    def test_DS (self):
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6)
        assert round(curve.DS(2e6), 9) == 100
        assert round(curve.DS(1.6e7), 9) == 50
        np.testing.assert_allclose(curve.DS(np.array([1.6e7, 2e6, 2.5e5])), [50, 100, 200])

    def test_damage (self):
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6)
        assert curve.damage(1e6, 100) == 0.5

    def test_cumulated_damage (self):
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6)
        mkvm = MarkovMatrix(range=np.array([100, 50]), mean=np.zeros(2), cycles=np.array([1e6, 1.6e6]))
        assert round(curve.cumulated_damage(mkvm), 12) == 0.6



class TestMultiSlopeFatigueCurve:

    def curve (self):
        return MultiSlopeFatigueCurve(
            SingleSlopeFatigueCurve(m=3, DS_ref=100, N_ref=2e6),
            SingleSlopeFatigueCurve(m=5, DS_ref=100, N_ref=2e6))

    def test_curves (self):
        curve = self.curve()
        assert len(curve.curves) == 2
        assert curve.curves[0].m == 3
        assert curve.curves[1].m == 5

    def test_N (self):
        curve = self.curve()
        assert curve.N(100) == 2e6
        assert curve.N(50) == 6.4e7     # governed by the m=5 slope
        assert curve.N(200) == 2.5e5    # governed by the m=3 slope
        np.testing.assert_allclose(curve.N(np.array([50., 100., 200.])), [6.4e7, 2e6, 2.5e5])

    def test_DS (self):
        curve = self.curve()
        assert round(curve.DS(2e6), 9) == 100
        assert round(curve.DS(2e7), 3) == round(100 * 0.1**(1/5), 3)    # governed by the m=5 slope
        assert round(curve.DS(2e5), 3) == round(100 * 10**(1/3), 3)     # governed by the m=3 slope
        np.testing.assert_allclose(curve.DS(np.array([6.4e7, 2e6, 2.5e5])), [50, 100, 200])



class TestDoubleSlopeFatigueCurve:

    def test_curves (self):
        curve = DoubleSlopeFatigueCurve(m1=3, m2=5, DS12=100, N12=2e6)
        assert isinstance(curve, MultiSlopeFatigueCurve)
        assert curve.curves == (
            SingleSlopeFatigueCurve(3, 100, 2e6),
            SingleSlopeFatigueCurve(5, 100, 2e6))

    def test_N (self):
        curve = DoubleSlopeFatigueCurve(m1=3, m2=5, DS12=100, N12=2e6)
        np.testing.assert_allclose(curve.N(np.array([50., 100., 200.])), [6.4e7, 2e6, 2.5e5])

    def test_DS (self):
        curve = DoubleSlopeFatigueCurve(m1=3, m2=5, DS12=100, N12=2e6)
        np.testing.assert_allclose(curve.DS(np.array([6.4e7, 2e6, 2.5e5])), [50, 100, 200])



class TestBoltFatigueCurve:

    def test_slopes (self):
        curve = BoltFatigueCurve(0.080)
        assert isinstance(curve, DoubleSlopeFatigueCurve)
        assert curve.curves[0].m == 3
        assert curve.curves[1].m == 5
        assert curve.curves[0].N_ref == 2e6
        assert curve.curves[1].N_ref == 2e6

    def test_size_effect (self):
        # No size effect up to 30 mm
        curve = BoltFatigueCurve(0.024)
        assert curve.curves[0].DS_ref == 50e6 / 1.1

        # First size effect factor between 30 mm and 72 mm
        curve = BoltFatigueCurve(0.048)
        assert curve.curves[0].DS_ref == pytest.approx(50e6 * (30/48)**0.1 / 1.1)

        # Both size effect factors above 72 mm
        curve = BoltFatigueCurve(0.080)
        assert curve.curves[0].DS_ref == pytest.approx(50e6 * (30/80)**0.1 * (72/80)**0.25 / 1.1)
        assert round(curve.curves[0].DS_ref / MPa, 2) == 40.14

    def test_custom_parameters (self):
        curve = BoltFatigueCurve(0.024, DS_ref=62*MPa, gamma_M=1.0)
        assert curve.curves[0].DS_ref == 62*MPa
        assert curve.curves[1].DS_ref == 62*MPa

    def test_N (self):
        curve = BoltFatigueCurve(0.024, DS_ref=50*MPa, gamma_M=1.0)
        assert curve.N(50*MPa) == pytest.approx(2e6)
        assert curve.N(100*MPa) == pytest.approx(2.5e5)
        assert curve.N(25*MPa) == pytest.approx(6.4e7)

    def test_DS (self):
        curve = BoltFatigueCurve(0.024, DS_ref=50*MPa, gamma_M=1.0)
        assert curve.DS(2e6) == pytest.approx(50*MPa)
        assert curve.DS(6.4e7) == pytest.approx(25*MPa)



class TestBoltFatigueAnalysis:

    def test_default_attributes (self):
        fseg = lflange_segment()
        mkvm = flange_markov_matrix()
        bfa = BoltFatigueAnalysis(fseg, mkvm)
        assert bfa.fseg is fseg
        assert bfa.flange_mkvm is mkvm
        assert bfa.custom_fatigue_curve is None
        assert bfa.allowable_damage == 1.0
        assert bfa.SMF == 1.0

    def test_fatigue_curve (self):
        fseg = lflange_segment()

        # Default curve
        bfa = BoltFatigueAnalysis(fseg, flange_markov_matrix())
        assert isinstance(bfa.fatigue_curve, BoltFatigueCurve)
        assert bfa.fatigue_curve.curves == BoltFatigueCurve(0.080).curves

        # Custom curve
        curve = SingleSlopeFatigueCurve(m=3, DS_ref=50*MPa, N_ref=2e6)
        bfa = BoltFatigueAnalysis(fseg, flange_markov_matrix(), custom_fatigue_curve=curve)
        assert bfa.fatigue_curve is curve

    def test_bolt_mkvm (self):
        fseg = lflange_segment()
        mkvm = flange_markov_matrix()
        bfa = BoltFatigueAnalysis(fseg, mkvm, SMF=1.2)

        bending_factor = 0.5 + 0.5*log(0.080/0.036) / log(150/36)
        expected = bolt_markov_matrix(fseg, mkvm, bending_factor, SMF=1.2)

        np.testing.assert_allclose(bfa.bolt_mkvm.range, expected.range)
        np.testing.assert_allclose(bfa.bolt_mkvm.mean, expected.mean)
        np.testing.assert_array_equal(bfa.bolt_mkvm.cycles, mkvm.cycles)
        assert bfa.bolt_mkvm.duration == mkvm.duration

    def test_damage (self):
        fseg = lflange_segment()
        bfa = BoltFatigueAnalysis(fseg, flange_markov_matrix())
        assert bfa.damage > 0
        assert bfa.damage == bfa.fatigue_curve.cumulated_damage(bfa.bolt_mkvm)

        # The damage grows with the stress multiplication factor
        assert BoltFatigueAnalysis(fseg, flange_markov_matrix(), SMF=1.2).damage > bfa.damage

    def test_fatigue_life (self):
        fseg = lflange_segment()
        bfa = BoltFatigueAnalysis(fseg, flange_markov_matrix(), allowable_damage=0.5)
        assert bfa.fatigue_life == pytest.approx(0.5 / bfa.damage * 20)
