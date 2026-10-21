import pytest
import os
import json
import copy
import warnings
from math import pi

import numpy as np
import openpyxl

import pyflange.io as io
from pyflange.bolts import StandardMetricBolt, ISOFlatWasher, ISOHexNut, RoundNut, MetricBolt, FlatWasher, HexNut, standard_hole_diameter
from pyflange.flangesegments import PolynomialLFlangeSegment, PolynomialLFlangeSegmentED2, Gap, shell_stiffness
from pyflange.fatigue import (MarkovMatrix, BoltFatigueCurve, SingleSlopeFatigueCurve,
                              DoubleSlopeFatigueCurve, BoltFatigueAnalysis)

deg = pi/180
MPa = 1e6

EXAMPLE_CASE_PATH = os.path.join(os.path.dirname(__file__), "..", "docsrc", "examples", "fatigue-case.json")


@pytest.fixture(autouse=True)
def random_seed ():
    # Make the random samples reproducible
    np.random.seed(20241021)


@pytest.fixture
def case ():
    # A complete case in the current format version (the documented example)
    with open(EXAMPLE_CASE_PATH) as file:
        return json.load(file)


@pytest.fixture
def legacy_case (case):
    # The same case in the unversioned 0.1.0 format, as written by Forrester
    legacy = copy.deepcopy(case)
    del legacy['format_version']
    del legacy['AllowableDamage']
    del legacy['SMF']
    del legacy['BoltFatigueCurve']['type']
    del legacy['Stats']['FatigueClassFactor']['DS_ref_mean']
    legacy['FlangeSegment']['bolt']['shank_diameter_ratio'] = 0.95125
    legacy['FlangeSegment']['k_shell'] = 2.5e9
    return legacy



class TestFormatVersion:

    def test_format_version (self, case, legacy_case):
        assert io.format_version(case) == "0.2.0"
        assert io.format_version(legacy_case) == "0.1.0"

    def test_upgrade_legacy_case (self, legacy_case):
        original = copy.deepcopy(legacy_case)
        upgraded = io.upgrade_case(legacy_case)
        assert upgraded['format_version'] == io.FORMAT_VERSION
        assert upgraded['Stats']['FatigueClassFactor']['DS_ref_mean'] == 62*MPa
        assert legacy_case == original     # input not modified

    def test_upgrade_legacy_case_without_stats (self, legacy_case):
        del legacy_case['Stats']
        assert 'Stats' not in io.upgrade_case(legacy_case)

    def test_upgrade_current_case (self, case):
        assert io.upgrade_case(case) == case

    def test_unsupported_version (self, case):
        case['format_version'] = "9.0.0"
        with pytest.raises(ValueError):
            io.upgrade_case(case)
        with pytest.raises(ValueError):
            io.create_fatiguecase(case)



class TestCreateBolt:

    def test_standard_metric_bolt (self):
        bolt = io.create_bolt(type="StandardMetricBolt", designation="M80", material_grade="10.9",
                              shank_length=0.16, shank_diameter_ratio=0.95, stud=True)
        assert bolt == StandardMetricBolt("M80", "10.9", shank_length=0.16, shank_diameter_ratio=0.95, stud=True)

    def test_pitch_shank_diameter_ratio (self):
        bolt = io.create_bolt(type="StandardMetricBolt", designation="M80", material_grade="10.9",
                              shank_length=0.16, shank_diameter_ratio="pitch", stud=True)
        assert bolt.shank_diameter == pytest.approx(bolt.thread_basic_pitch_diameter)

    def test_metric_bolt (self):
        params = dict(nominal_diameter=0.08, thread_pitch=0.006, yield_stress=900e6, ultimate_tensile_stress=1000e6,
                      elastic_modulus=200e9, poissons_ratio=0.29, shank_length=0.2, shank_diameter_ratio=0.9, stud=False)
        assert io.create_bolt(type="MetricBolt", **params) == MetricBolt(**params)

    def test_unsupported_type (self):
        with pytest.raises(ValueError):
            io.create_bolt(type="ImperialBolt")



class TestCreateWasher:

    def test_no_washer (self):
        assert io.create_washer(None) is None

    def test_iso_flat_washer (self):
        assert io.create_washer({"type": "ISOFlatWasher", "designation": "M64"}) == ISOFlatWasher("M64")

    def test_flat_washer (self):
        params = dict(outer_diameter=0.12, inner_diameter=0.07, thickness=0.01, elastic_modulus=200e9, poissons_ratio=0.3)
        assert io.create_washer({"type": "FlatWasher", **params}) == FlatWasher(**params)

    def test_unsupported_type (self):
        with pytest.raises(ValueError):
            io.create_washer({"type": "SpringWasher"})



class TestCreateNut:

    def test_standard_nuts (self):
        assert io.create_nut(type="ISOHexNut", designation="M64") == ISOHexNut("M64")
        assert io.create_nut(type="RoundNut", designation="M80") == RoundNut("M80")

    def test_hex_nut (self):
        params = dict(nominal_diameter=0.08, thickness=0.064, inscribed_diameter=0.115, circumscribed_diameter=0.1275,
                      bearing_diameter=0.140, elastic_modulus=200e9, poissons_ratio=0.3)
        assert io.create_nut(type="HexNut", **params) == HexNut(**params)

    def test_unsupported_type (self):
        with pytest.raises(ValueError):
            io.create_nut(type="WingNut")



class TestCreateFlangeSegment:

    def test_example_case (self, case):
        params = case['FlangeSegment']
        fseg = io.create_flangesegment(**params)
        assert type(fseg) is PolynomialLFlangeSegment
        assert fseg.a == 0.150
        assert fseg.bolt == StandardMetricBolt("M80", "10.9", shank_length=0.16, shank_diameter_ratio="pitch", stud=True)
        assert fseg.washer is None
        assert fseg.nut == RoundNut("M80")
        assert fseg.gap == Gap(0.0013, 30*deg, 1.0)
        assert fseg.Do == standard_hole_diameter("M80")
        assert fseg.Fv == pytest.approx(0.75301 * fseg.bolt.thread_cross_section.area * fseg.bolt.yield_stress)

    def test_default_type (self, case):
        del case['FlangeSegment']['type']
        assert io.DEFAULT_FLANGE_SEGMENT_TYPE == "PolynomialLFlangeSegment"
        assert type(io.create_flangesegment(**case['FlangeSegment'])) is PolynomialLFlangeSegment

    def test_ed2_type (self, case):
        case['FlangeSegment']['type'] = "PolynomialLFlangeSegmentED2"
        assert type(io.create_flangesegment(**case['FlangeSegment'])) is PolynomialLFlangeSegmentED2

    def test_unsupported_type (self, case):
        case['FlangeSegment']['type'] = "PolynomialTFlangeSegment"
        with pytest.raises(ValueError):
            io.create_flangesegment(**case['FlangeSegment'])

    def test_absolute_preload (self, case):
        params = case['FlangeSegment']
        fseg1 = io.create_flangesegment(**params)
        del params['preload_ratio']
        params['Fv'] = fseg1.Fv
        fseg2 = io.create_flangesegment(**params)
        assert fseg2.Fv == fseg1.Fv
        assert fseg2.bolt_preload == fseg1.bolt_preload

    def test_preload_specified_once (self, case):
        params = case['FlangeSegment']
        params['Fv'] = 2.8e6
        with pytest.raises(ValueError):
            io.create_flangesegment(**params)        # both preload_ratio and Fv
        del params['Fv']
        del params['preload_ratio']
        with pytest.raises(ValueError):
            io.create_flangesegment(**params)        # neither

    def test_shell_stiffness (self, case):
        params = case['FlangeSegment']
        fseg = io.create_flangesegment(**params)
        assert fseg.k_shell == pytest.approx(float(shell_stiffness(params['R'], params['s'], params['gap']['angle'])))

        params['k_shell'] = 3e9
        assert io.create_flangesegment(**params).k_shell == 3e9

        params['k_shell'] = None
        assert io.create_flangesegment(**params).k_shell is None

    def test_extra_keys_ignored (self, case):
        params = case['FlangeSegment']
        params['BoltPreload'] = 1.0
        params['Point1'] = {"Z": 0}
        assert isinstance(io.create_flangesegment(**params), PolynomialLFlangeSegment)



class TestCreateFatigueCurve:

    def test_bolt_fatigue_curve (self):
        curve = io.create_fatiguecurve(BoltNominalDiameter=0.08, DS_ref=50*MPa, gamma_M=1.1)
        assert isinstance(curve, BoltFatigueCurve)
        assert curve.N(40*MPa) == BoltFatigueCurve(0.08, 50*MPa, 1.1).N(40*MPa)

        curve = io.create_fatiguecurve(type="BoltFatigueCurve", BoltNominalDiameter=0.08, DS_ref=50*MPa, gamma_M=1.1)
        assert isinstance(curve, BoltFatigueCurve)

    def test_other_curves (self):
        curve = io.create_fatiguecurve(type="SingleSlopeFatigueCurve", m=3, DS_ref=100*MPa, N_ref=2e6)
        assert curve == SingleSlopeFatigueCurve(3, 100*MPa, 2e6)

        curve = io.create_fatiguecurve(type="DoubleSlopeFatigueCurve", m1=3, m2=5, DS_ref=40*MPa, N_ref=2e6)
        assert isinstance(curve, DoubleSlopeFatigueCurve)
        assert curve.N(40*MPa) == pytest.approx(2e6)

    def test_unsupported_type (self):
        with pytest.raises(ValueError):
            io.create_fatiguecurve(type="TripleSlopeFatigueCurve")



def test_create_markovmatrix ():
    mkvm = io.create_markovmatrix(range=[1, 2], mean=[3, 4], cycles=[5, 6], duration=20)
    assert isinstance(mkvm, MarkovMatrix)
    np.testing.assert_array_equal(mkvm.range, [1., 2.])
    np.testing.assert_array_equal(mkvm.mean, [3., 4.])
    np.testing.assert_array_equal(mkvm.cycles, [5., 6.])
    assert mkvm.duration == 20



class TestFatigueCase:

    def test_load_fatiguecase (self):
        fcase = io.load_fatiguecase(EXAMPLE_CASE_PATH)
        assert isinstance(fcase, BoltFatigueAnalysis)
        assert isinstance(fcase.fatigue_curve, BoltFatigueCurve)
        assert fcase.allowable_damage == 1.0
        assert fcase.SMF == 1.0
        assert 0 < fcase.damage < np.inf

    def test_smf (self, case):
        case['SMF'] = 1.2
        assert io.create_fatiguecase(case).SMF == 1.2
        assert io.create_fatiguecase(case, SMF=1.5).SMF == 1.5
        del case['SMF']
        assert io.create_fatiguecase(case).SMF == 1.0

    def test_allowable_damage (self, case):
        case['AllowableDamage'] = 0.5
        assert io.create_fatiguecase(case).allowable_damage == 0.5
        del case['AllowableDamage']
        assert io.create_fatiguecase(case).allowable_damage == 1.0

    def test_legacy_case (self, case, legacy_case):
        # A 0.1.0 case produces the same results as the equivalent 0.2.0 case
        case['FlangeSegment']['bolt']['shank_diameter_ratio'] = 0.95125
        case['FlangeSegment']['k_shell'] = 2.5e9
        assert io.create_fatiguecase(legacy_case).damage == io.create_fatiguecase(case).damage

    def test_no_warnings (self, case):
        case['MarkovMatrix']['range'][0] = 0.0     # zero-range bin
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert io.create_fatiguecase(case).damage > 0



class TestDistributions:

    def test_normal (self):
        dist = io.create_distribution(Distribution="Normal", Mean=10.0, CoV=0.2)
        assert dist.mean() == pytest.approx(10.0)
        assert dist.std() == pytest.approx(2.0)

    def test_lognormal (self):
        dist = io.create_distribution(Distribution="LogNormal", Mean=10.0, CoV=0.2)
        assert dist.mean() == pytest.approx(10.0)
        assert dist.std() == pytest.approx(2.0)

    def test_weibull (self):
        dist = io.create_distribution(Distribution="Weibull", Shape=1.0, Scale=2.0)
        assert dist.mean() == pytest.approx(2.0)

    def test_beta (self):
        dist = io.create_distribution(Distribution="Beta", Alpha=2, Beta=2)
        assert dist.support() == (0, 1)
        dist = io.create_distribution(Distribution="Beta", Alpha=2, Beta=2, x0=1, x1=3)
        assert dist.support() == (1, 3)
        assert dist.mean() == pytest.approx(2.0)

    def test_unsupported (self):
        with pytest.raises(ValueError):
            io.create_distribution(Distribution="Cauchy")
        with pytest.raises(ValueError):
            io.create_sampler(Distribution="Cauchy")

    def test_samplers (self):
        values = [next(io.create_sampler(Distribution="Normal", Mean=10.0, CoV=0.0)) for i in range(3)]
        assert values == [10.0, 10.0, 10.0]

        samp = io.create_sampler(Distribution="Product",
                                 Dist1={"Distribution": "Normal", "Mean": 2.0, "CoV": 0.0},
                                 Dist2={"Distribution": "Normal", "Mean": 3.0, "CoV": 0.0})
        assert next(samp) == 6.0



class TestGapSampler:

    def test_amd1 (self, case):
        samp = io.create_gap_sampler(8.0, **case['Stats']['Gap'])
        gap = next(samp)
        assert isinstance(gap, Gap)
        assert 0 < gap.angle < pi

    def test_ed2 (self):
        samp = io.create_gap_sampler(8.0, JointDistribution="IEC 61400-6:2020 ED2", FlatnessTolerance=0.0014)
        gap = next(samp)
        assert isinstance(gap, Gap)
        assert gap.shape_factor == 1.0

    def test_clayton (self):
        samp = io.create_gap_sampler(8.0, JointDistribution="ClaytonCopula", Theta=2.0,
                    Height={"Distribution": "LogNormal", "Mean": 0.001, "CoV": 0.8},
                    Angle={"Distribution": "Beta", "Alpha": 2, "Beta": 3, "x0": 0, "x1": pi},
                    ShapeFactor={"Distribution": "Normal", "Mean": 1.0, "CoV": 0.15})
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            gap = next(samp)
        assert isinstance(gap, Gap)

    def test_unsupported (self):
        with pytest.raises(ValueError):
            io.create_gap_sampler(8.0, JointDistribution="Uniform")



class TestFatigueCaseSampler:

    def test_load_fatiguecase_sampler (self):
        samp = io.load_fatiguecase_sampler(EXAMPLE_CASE_PATH)
        fcase = next(samp)
        assert isinstance(fcase, BoltFatigueAnalysis)
        assert type(fcase.fseg) is PolynomialLFlangeSegment
        assert fcase.fseg.PPL == 0 and fcase.fseg.DFT == 0
        assert 0 < fcase.damage < np.inf

    def test_segment_type (self, case):
        # ED2 cases generate ED2 flange segments
        case['FlangeSegment']['type'] = "PolynomialLFlangeSegmentED2"
        case['Stats']['Gap'] = {"JointDistribution": "IEC 61400-6:2020 ED2", "FlatnessTolerance": 0.0014}
        data = io.upgrade_case(case)
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        samp = io.create_fatiguecase_sampler(fseg, io.create_markovmatrix(**data['MarkovMatrix']), **data['Stats'])
        assert type(next(samp).fseg) is PolynomialLFlangeSegmentED2

        # Legacy (0.1.0) ED2 cases generate ED2 flange segments too
        del case['format_version']
        del case['Stats']['FatigueClassFactor']['DS_ref_mean']
        data = io.upgrade_case(case)
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        samp = io.create_fatiguecase_sampler(fseg, io.create_markovmatrix(**data['MarkovMatrix']), **data['Stats'])
        assert type(next(samp).fseg) is PolynomialLFlangeSegmentED2

    def test_default_segment_type (self, case):
        # Without type, AMD1 flange segments are generated
        del case['FlangeSegment']['type']
        assert next(io.create_fatiguecase_sampler(
            io.create_flangesegment(**case['FlangeSegment']),
            io.create_markovmatrix(**case['MarkovMatrix']),
            **case['Stats'])).fseg.__class__ is PolynomialLFlangeSegment

    def test_absolute_mean_preload (self, case):
        data = io.upgrade_case(case)
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        mkvm = io.create_markovmatrix(**data['MarkovMatrix'])
        data['Stats']['Preload'] = {"Distribution": "Normal", "Mean": 2.9e6, "CoV": 0.0}
        data['Stats']['PlasticLoss'] = {"Distribution": "Normal", "Mean": 0.0, "CoV": 0.0}
        data['Stats']['DFT'] = {"Distribution": "Normal", "Mean": 0.0, "CoV": 0.0}
        fcase = next(io.create_fatiguecase_sampler(fseg, mkvm, **data['Stats']))
        assert fcase.fseg.Fv == pytest.approx(2.9e6)

    def test_constant_smf (self, case):
        data = io.upgrade_case(case)
        data['Stats']['StressMultiplicationFactor'] = 1.3
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        mkvm = io.create_markovmatrix(**data['MarkovMatrix'])
        assert next(io.create_fatiguecase_sampler(fseg, mkvm, **data['Stats'])).SMF == 1.3

    def test_missing_DS_ref_mean (self, case):
        data = io.upgrade_case(case)
        del data['Stats']['FatigueClassFactor']['DS_ref_mean']
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        mkvm = io.create_markovmatrix(**data['MarkovMatrix'])
        with pytest.raises(KeyError):
            io.create_fatiguecase_sampler(fseg, mkvm, **data['Stats'])

    def test_fatigue_class (self, case):
        data = io.upgrade_case(case)
        data['Stats']['FatigueClassFactor'] = {"Distribution": "Normal", "Mean": 1.0, "CoV": 0.0, "DS_ref_mean": 70*MPa}
        fseg = io.create_flangesegment(**data['FlangeSegment'])
        mkvm = io.create_markovmatrix(**data['MarkovMatrix'])
        fcase = next(io.create_fatiguecase_sampler(fseg, mkvm, **data['Stats']))
        assert fcase.fatigue_curve.DS_ref == pytest.approx(70*MPa)
        assert fcase.fatigue_curve.gamma_M == 1.0



class TestMonteCarlo:

    def test_sample_fatiguecase (self, case):
        progress = []
        sample = io.sample_fatiguecase(case, 5, seed=1, progress=progress.append)
        assert list(sample.columns) == io.SAMPLE_COLUMNS
        assert len(sample) == 5
        assert progress == [1, 2, 3, 4, 5]
        assert np.all(sample.Utilization == sample.Damage / sample.AllowableDamage)

    def test_seed (self, case):
        sample1 = io.sample_fatiguecase(case, 3, seed=1)
        sample2 = io.sample_fatiguecase(case, 3, seed=1)
        sample3 = io.sample_fatiguecase(case, 3, seed=2)
        assert sample1.equals(sample2)
        assert not sample1.equals(sample3)

    def test_legacy_case (self, case, legacy_case):
        # Same seed, same sample, whether the case is in 0.1.0 or 0.2.0 format
        case['FlangeSegment']['bolt']['shank_diameter_ratio'] = 0.95125
        assert io.sample_fatiguecase(legacy_case, 3, seed=1).equals(io.sample_fatiguecase(case, 3, seed=1))

    def test_moment_range_coefficient (self, case):
        # The coefficient is evaluated on the first non-zero range
        case['MarkovMatrix']['range'][0] = 0.0
        sample = io.sample_fatiguecase(case, 3, seed=1)
        assert np.all(np.isfinite(sample.MomentRangeCoefficient))

    @pytest.mark.parametrize("model", ["AMD1", "ED2", "Clayton"])
    def test_create_realization (self, case, model):
        if model == "ED2":
            case['FlangeSegment']['type'] = "PolynomialLFlangeSegmentED2"
            case['Stats']['Gap'] = {"JointDistribution": "IEC 61400-6:2020 ED2", "FlatnessTolerance": 0.0014}
        elif model == "Clayton":
            case['Stats']['Gap'] = {"JointDistribution": "ClaytonCopula", "Theta": 2.0,
                    "Height": {"Distribution": "LogNormal", "Mean": 0.001, "CoV": 0.8},
                    "Angle": {"Distribution": "Beta", "Alpha": 2, "Beta": 3, "x0": 0.2, "x1": pi},
                    "ShapeFactor": {"Distribution": "Normal", "Mean": 1.0, "CoV": 0.15}}

        # Each rebuilt realization gives the same results of the sampled one
        sample = io.sample_fatiguecase(case, 4, seed=3)
        for i in range(len(sample)):
            fcase = io.create_realization(case, sample.loc[i])
            assert fcase.fseg.bolt_preload == pytest.approx(sample.BoltPreload[i], rel=1e-12)
            assert fcase.damage == pytest.approx(sample.Damage[i], rel=1e-9)
            assert fcase.fatigue_life == pytest.approx(sample.FatigueLife[i], rel=1e-9)



class TestExport:

    def test_bolt_to_dict (self):
        bolt = StandardMetricBolt("M80", "10.9", shank_length=0.16, shank_diameter_ratio=0.95, stud=True)
        d = io.bolt_to_dict(bolt)
        assert d['type'] == "MetricBolt"
        assert io.create_bolt(**d) == bolt
        with pytest.raises(ValueError):
            io.bolt_to_dict("M80")

    def test_washer_to_dict (self):
        assert io.washer_to_dict(None) is None
        washer = ISOFlatWasher("M64")
        assert io.create_washer(io.washer_to_dict(washer)) == washer
        with pytest.raises(ValueError):
            io.washer_to_dict("M64")

    def test_nut_to_dict (self):
        nut = RoundNut("M80")
        assert io.create_nut(**io.nut_to_dict(nut)) == nut
        with pytest.raises(ValueError):
            io.nut_to_dict("M80")

    def test_gap_to_dict (self):
        assert io.gap_to_dict(Gap(0.001, 0.5, 1.1)) == {"height": 0.001, "angle": 0.5, "shape_factor": 1.1}
        with pytest.raises(ValueError):
            io.gap_to_dict(None)

    def test_fatiguecurve_to_dict (self):
        curve = BoltFatigueCurve(0.08, 62*MPa, 1.0)
        d = io.fatiguecurve_to_dict(curve)
        assert d['type'] == "BoltFatigueCurve"
        assert d['DS12'] == curve.curves[0].DS_ref
        assert io.create_fatiguecurve(**d).N(40*MPa) == curve.N(40*MPa)

        curve = SingleSlopeFatigueCurve(3, 100*MPa, 2e6)
        assert io.create_fatiguecurve(**io.fatiguecurve_to_dict(curve)) == curve

        curve = DoubleSlopeFatigueCurve(3, 5, 40*MPa, 2e6)
        assert io.create_fatiguecurve(**io.fatiguecurve_to_dict(curve)).N(30*MPa) == curve.N(30*MPa)

        with pytest.raises(ValueError):
            io.fatiguecurve_to_dict(None)

    def test_markovmatrix_to_dict (self):
        with pytest.raises(ValueError):
            io.markovmatrix_to_dict(None)

    def test_flangesegment_to_dict (self, case):
        case['FlangeSegment']['type'] = "PolynomialLFlangeSegmentED2"
        fseg = io.create_flangesegment(**case['FlangeSegment'])
        d = io.flangesegment_to_dict(fseg)
        assert d['type'] == "PolynomialLFlangeSegmentED2"
        assert 'Point1' not in d

        d = io.flangesegment_to_dict(fseg, verbose=True)
        assert d['BoltPreload'] == pytest.approx(fseg.bolt_preload)
        assert len(d['TransferFunctions']['Z']) == 100

        with pytest.raises(ValueError):
            io.flangesegment_to_dict(None)

    @pytest.mark.parametrize("verbose", [False, True])
    def test_fatiguecase_round_trip (self, case, verbose):
        fcase = io.create_fatiguecase(case)
        d = io.fatiguecase_to_dict(fcase, verbose=verbose)
        assert d['format_version'] == io.FORMAT_VERSION
        assert d['damage'] == fcase.damage
        assert ('FatigueLife' in d) == verbose

        # The output is JSON-serializable and can be loaded back
        d = json.loads(json.dumps(d))
        assert io.create_fatiguecase(d).damage == pytest.approx(fcase.damage, rel=1e-12)

        with pytest.raises(ValueError):
            io.fatiguecase_to_dict(None)

    def test_dict_to_json (self, tmp_path):
        path = tmp_path / "data.json"
        io.dict_to_json({"a": np.float64(1.5), "b": np.array([1, 2]), "c": {"d": np.bool_(True)}}, path)
        with open(path) as file:
            assert json.load(file) == {"a": 1.5, "b": [1, 2], "c": {"d": True}}

    def test_dict_to_excel (self, tmp_path):
        path = tmp_path / "data.xlsx"
        io.dict_to_excel({"a": 10, "b": {"c": np.float64(20.5), "d": [1, 2]}, "e": None}, path)
        ws = openpyxl.load_workbook(path)["CaseData"]
        assert [[cell.value for cell in row] for row in ws.iter_rows()] == [
            ["KEY", "VALUE"], ["a", 10], ["b.c", 20.5], ["b.d", "1, 2"], ["e", None]]
        assert ws.defined_names["b.c"].attr_text == "CaseData!$B$3"
