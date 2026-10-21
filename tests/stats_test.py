import pytest
import pyflange.stats as stats
from pyflange.flangesegments import Gap, PolynomialLFlangeSegment, shell_stiffness
from pyflange.fatigue import MarkovMatrix, BoltFatigueCurve, BoltFatigueAnalysis
from pyflange.bolts import MetricBolt, HexNut
from math import pi, gamma
from itertools import repeat
from types import GeneratorType
import scipy.stats
import numpy as np

deg = pi/180
MPa = 1e6


@pytest.fixture(autouse=True)
def random_seed ():
    # Make the random samples reproducible
    np.random.seed(20241021)


def take (samp, n):
    return np.array([next(samp) for i in range(n)])


def lflange_segment_params ():
    return dict(
        a = 0.2325,
        b = 0.1665,
        s = 0.0720,
        t = 0.2000,
        R = 3.75,
        central_angle = 2*pi/120,
        Zg = -14795000/120,
        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            shank_diameter_ratio = 76.1/80,
            shank_length = 0.270,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            stud = True),
        Do = 0.086,
        washer = None,
        nut = HexNut(
            nominal_diameter = 0.080,
            thickness = 0.064,
            inscribed_diameter = 0.115,
            circumscribed_diameter = 0.1275,
            bearing_diameter = 0.140),
        s_ratio = 100/72)



def test_gap_height_distribution ():
    D = 7.5
    u = 0.0014
    L = pi/6 * D/2

    gap_dist = stats.gap_height_distribution(D, u, L)

    assert round(gap_dist.mean(), 6) == 0.000288
    assert round(gap_dist.std(), 6) == 0.000350
    assert round(gap_dist.ppf(0.95), 6) == 0.000876



def test_sampler ():
    samp = stats.sampler(scipy.stats.uniform(2, 3))
    assert isinstance(samp, GeneratorType)

    values = take(samp, 5000)
    assert np.all((values >= 2) & (values <= 5))
    assert np.mean(values) == pytest.approx(3.5, rel=0.02)



def test_norm_sampler ():
    values = take(stats.norm_sampler(12.0, 0.25), 10000)
    assert np.mean(values) == pytest.approx(12.0, rel=0.02)
    assert np.std(values) == pytest.approx(3.0, rel=0.05)



def test_lognorm_sampler ():
    values = take(stats.lognorm_sampler(12.0, 0.25), 10000)
    assert np.all(values > 0)
    assert np.mean(values) == pytest.approx(12.0, rel=0.02)
    assert np.std(values) == pytest.approx(3.0, rel=0.05)



def test_weibull_sampler ():
    values = take(stats.weibull_sampler(1.5, 2.0), 10000)
    assert np.all(values > 0)
    assert np.mean(values) == pytest.approx(2.0 * gamma(1 + 1/1.5), rel=0.03)

    # Default scale
    values = take(stats.weibull_sampler(1.5), 10000)
    assert np.mean(values) == pytest.approx(gamma(1 + 1/1.5), rel=0.03)



def test_beta_sampler ():
    values = take(stats.beta_sampler(2, 5, 1, 3), 10000)
    assert np.all((values >= 1) & (values <= 3))
    assert np.mean(values) == pytest.approx(1 + 2 * 2/7, rel=0.02)

    # Default interval [0, 1]
    values = take(stats.beta_sampler(2, 5), 10000)
    assert np.all((values >= 0) & (values <= 1))
    assert np.mean(values) == pytest.approx(2/7, rel=0.03)



def test_standard_gap_sampler ():
    D = 7.5
    u = 0.0014

    # Default samplers
    gaps = take(stats.standard_gap_sampler(D, u), 200)
    for gap in gaps:
        assert isinstance(gap, Gap)
        assert 0 < gap.angle < pi
        assert gap.height <= 2.0 * stats.gap_height_distribution(D, u, gap.angle*D/2).ppf(0.95)

    # Custom samplers and gap height cap
    p95 = stats.gap_height_distribution(D, u, 30*deg*D/2).ppf(0.95)
    gaps = take(stats.standard_gap_sampler(D, u,
                    gap_angle_sampler = repeat(30*deg),
                    gap_shape_factor_sampler = repeat(1.2),
                    gap_height_cap = 1.0), 200)
    for gap in gaps:
        assert gap.angle == 30*deg
        assert gap.shape_factor == 1.2
        assert 0 < gap.height <= p95

    # Gap angles outside (0, pi) are discarded
    angles = iter([200*deg, 60*deg])
    gap = next(stats.standard_gap_sampler(D, u, gap_angle_sampler=angles, gap_height_cap=1e9))
    assert gap.angle == 60*deg



def test_standard_PolynomialLFlangeSegment_sampler ():
    params = lflange_segment_params()
    gap = Gap(height=0.0005, angle=30*deg, shape_factor=1.0)

    samp = stats.standard_PolynomialLFlangeSegment_sampler(
        **params,
        preload_sampler = repeat(2800e3),
        ppl_sampler = repeat(0.05),
        dft_sampler = repeat(300e-6),
        gap_sampler = repeat(gap),
        tilt_sampler = repeat(-0.5*deg),
        settlement_factor = 0.5)

    fseg = next(samp)
    assert isinstance(fseg, PolynomialLFlangeSegment)
    assert fseg.gap is gap
    assert fseg.tilt_angle == pytest.approx(0.5*deg)    # reduced to [0, pi]
    assert fseg.k_shell == pytest.approx(shell_stiffness(params['R'], params['s'], gap.angle), rel=1e-2)

    # The preload losses are already included in the sampled preload
    assert fseg.PPL == 0
    assert fseg.DFT == 0
    reference_fseg = PolynomialLFlangeSegment(**params, gap=gap, Fv=2800e3, PPL=0.05, DFT=300e-6,
                                              settlement_factor=0.5, k_shell=fseg.k_shell)
    assert fseg.Fv == pytest.approx(reference_fseg.bolt_preload)
    assert fseg.bolt_preload == pytest.approx(reference_fseg.bolt_preload)

    # Custom shell stiffness
    samp = stats.standard_PolynomialLFlangeSegment_sampler(
        **params,
        preload_sampler = repeat(2800e3),
        ppl_sampler = repeat(0.0),
        dft_sampler = repeat(0.0),
        gap_sampler = repeat(gap),
        k_shell = 1.5e9)
    fseg = next(samp)
    assert fseg.k_shell == 1.5e9
    assert 0 <= fseg.tilt_angle <= pi

    # The preload is averaged over half of the bolts in the gap
    gap = Gap(height=0.0005, angle=60*deg)
    preloads = iter(np.arange(1, 100) * 1e5)
    samp = stats.standard_PolynomialLFlangeSegment_sampler(
        **params,
        preload_sampler = preloads,
        ppl_sampler = repeat(0.0),
        dft_sampler = repeat(0.0),
        gap_sampler = repeat(gap))
    fseg = next(samp)
    n = round(gap.angle / params['central_angle'] / 2)
    assert n == 10
    assert fseg.Fv == pytest.approx(np.mean(np.arange(1, n+1) * 1e5))



def test_standard_markov_matrix_sampler ():
    mkvm = MarkovMatrix(
        range = np.array([1.0, 2.0, 3.0]),
        mean = np.array([4.0, 5.0, 6.0]),
        cycles = np.array([10, 20, 30]),
        duration = 20)

    samp = stats.standard_markov_matrix_sampler(mkvm, mean_range_coeff=1.1, range_CoV=0.12)
    coeffs = []
    for i in range(5000):
        rnd_mkvm = next(samp)
        assert isinstance(rnd_mkvm, MarkovMatrix)
        np.testing.assert_array_equal(rnd_mkvm.cycles, mkvm.cycles)
        np.testing.assert_array_equal(rnd_mkvm.mean, mkvm.mean)
        assert rnd_mkvm.duration == mkvm.duration

        # All the ranges are scaled by the same coefficient
        ratios = rnd_mkvm.range / mkvm.range
        np.testing.assert_allclose(ratios, ratios[0])
        coeffs.append(ratios[0])

    assert np.mean(coeffs) == pytest.approx(1.1, rel=0.02)
    assert np.std(coeffs) / np.mean(coeffs) == pytest.approx(0.12, rel=0.05)



def test_standard_bolt_fatigue_curve_sampler ():
    d = 0.080
    reference_DS = BoltFatigueCurve(d, 62*MPa, gamma_M=1.0).curves[0].DS_ref

    samp = stats.standard_bolt_fatigue_curve_sampler(d, mean_stress_factor=1.05, stress_factor_CoV=0.10)
    coeffs = []
    for i in range(5000):
        curve = next(samp)
        assert isinstance(curve, BoltFatigueCurve)
        assert curve.curves[0].DS_ref == curve.curves[1].DS_ref
        coeffs.append(curve.curves[0].DS_ref / reference_DS)

    assert np.mean(coeffs) == pytest.approx(1.05, rel=0.02)
    assert np.std(coeffs) / np.mean(coeffs) == pytest.approx(0.10, rel=0.05)

    # Custom mean reference stress range
    samp = stats.standard_bolt_fatigue_curve_sampler(0.024, stress_factor_CoV=1e-9, DS_ref_mean=50*MPa)
    assert next(samp).curves[0].DS_ref == pytest.approx(50*MPa)



def test_fatigue_case_sampler ():
    params = lflange_segment_params()
    fseg = PolynomialLFlangeSegment(**params, gap=Gap(0.0005, 30*deg), Fv=2800e3)
    mkvm = MarkovMatrix(range=np.array([1e6]), mean=np.array([1e7]), cycles=np.array([1e6]))
    curve = BoltFatigueCurve(0.080)

    # Constant SMF
    samp = stats.fatigue_case_sampler(repeat(fseg), repeat(mkvm), repeat(curve), repeat(0.8), SMF=1.2)
    fcase = next(samp)
    assert isinstance(fcase, BoltFatigueAnalysis)
    assert fcase.fseg is fseg
    assert fcase.flange_mkvm is mkvm
    assert fcase.custom_fatigue_curve is curve
    assert fcase.allowable_damage == 0.8
    assert fcase.SMF == 1.2

    # Default SMF
    samp = stats.fatigue_case_sampler(repeat(fseg), repeat(mkvm), repeat(curve), repeat(1.0))
    assert next(samp).SMF == 1.0

    # Random SMF
    def SMF_sampler ():
        yield 1.1
        yield 1.3
    samp = stats.fatigue_case_sampler(repeat(fseg), repeat(mkvm), repeat(curve), repeat(1.0), SMF=SMF_sampler())
    assert next(samp).SMF == 1.1
    assert next(samp).SMF == 1.3
