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

"""Reading and writing flange fatigue cases in the PyFlange JSON format.

This module implements the PyFlange JSON case format, documented in detail in
the [JSON format reference](../json-format.md). In particular it provides:

- Loaders, that convert a JSON case file (or the equivalent dictionary) into
  PyFlange objects: `load_fatiguecase`, `load_fatiguecase_sampler` and the
  lower-level `create_*` functions.
- Exporters, that convert PyFlange objects back into dictionaries in the same
  format (`*_to_dict` functions) and write dictionaries to JSON or Excel files
  (`dict_to_json`, `dict_to_excel`).
- Monte Carlo helpers: `sample_fatiguecase`, that runs a simulation and returns
  a table of realizations, and `create_realization`, that rebuilds a single
  realization from that table.

Format versions:
    The current format version is `0.2.0` (`FORMAT_VERSION`). Files without a
    `format_version` field are assumed to be version `0.1.0` (the format used
    by the former Forrester tool) and are upgraded on the fly by
    `upgrade_case`, so they keep producing the same results.

Units:
    All values are expressed in SI units (m, N, Pa, rad), except where stated
    otherwise in the format reference.

Example:
    ```python
    from pyflange.io import load_fatiguecase, fatiguecase_to_dict, dict_to_json

    fcase = load_fatiguecase("case.json")       # a BoltFatigueAnalysis object
    print(fcase.damage, fcase.fatigue_life)

    dict_to_json(fatiguecase_to_dict(fcase, verbose=True), "case-res.json")
    ```
"""

import copy
import json

import numpy as np


FORMAT_VERSION = "0.2.0"
"""str: The version of the JSON format written by this module."""

SUPPORTED_FORMAT_VERSIONS = ("0.1.0", "0.2.0")
"""tuple: The JSON format versions that this module can read."""

LEGACY_DS_REF_MEAN = 62e6
"""float: Mean fatigue reference stress range implicitly used by format 0.1.0
samplers, in Pa. Format 0.1.0 files get it assigned when upgraded."""

FLANGE_SEGMENT_TYPES = ("PolynomialLFlangeSegment", "PolynomialLFlangeSegmentED2")
"""tuple: The flange segment types supported by the JSON format."""

DEFAULT_FLANGE_SEGMENT_TYPE = "PolynomialLFlangeSegment"
"""str: The flange segment type used when `FlangeSegment.type` is omitted
(IEC 61400-6 AMD1 model)."""

SAMPLE_COLUMNS = ['BoltPreload', 'GapShapeFactor', 'TiltAngle', 'GapAngle', 'GapHeight',
                  'FatigueClass', 'MomentRangeCoefficient', 'Damage', 'AllowableDamage',
                  'FatigueLife', 'SMF', 'Utilization']
"""list: The columns of the table returned by `sample_fatiguecase`."""





# -----------------------------------------------------------------------------
#   FORMAT VERSIONS
# -----------------------------------------------------------------------------

def format_version(data):
    """Returns the format version of a case dictionary.

    Args:
        data (dict): A case dictionary, as read from a JSON case file.

    Returns:
        str: The value of the `format_version` field, or `"0.1.0"` if the
            field is missing.
    """
    return str(data.get("format_version", "0.1.0"))


def upgrade_case(data):
    """Converts a case dictionary to the current format version.

    The input dictionary is not modified: a converted copy is returned.
    Upgrading a version `0.1.0` case:

    - sets `format_version` to the current version;
    - adds `DS_ref_mean = 62 MPa` to `Stats.FatigueClassFactor`, if the
      `Stats` section is present, which is the value implicitly used by the
      version `0.1.0` Monte Carlo sampler.

    Args:
        data (dict): A case dictionary in any supported format version.

    Returns:
        dict: An equivalent case dictionary in the current format version.

    Raises:
        ValueError: If the format version is not supported.
    """
    version = format_version(data)
    if version not in SUPPORTED_FORMAT_VERSIONS:
        raise ValueError(f"Unsupported case format version '{version}'. "
                         f"Supported versions: {', '.join(SUPPORTED_FORMAT_VERSIONS)}.")

    data = copy.deepcopy(data)

    if version == "0.1.0" and "Stats" in data:
        fcf = data["Stats"].get("FatigueClassFactor")
        if fcf is not None and "DS_ref_mean" not in fcf:
            fcf["DS_ref_mean"] = LEGACY_DS_REF_MEAN

    data["format_version"] = FORMAT_VERSION
    return data





# -----------------------------------------------------------------------------
#   IMPORT DATA
# -----------------------------------------------------------------------------

def load_inputfile(file_path):
    """Loads a JSON input file to a dictionary.

    The content is returned as it is, without any format upgrade or
    validation. Use `upgrade_case` to convert it to the current version.

    Args:
        file_path (str): Path of the JSON file.

    Returns:
        dict: The content of the JSON file.
    """
    with open(file_path, "r") as file:
        return json.load(file)


def create_bolt(**params):
    """Creates a bolt object from a dictionary of parameters.

    Args:
        **params (Any): The bolt parameters. The `type` key selects the bolt
            kind: `"StandardMetricBolt"` (keys `designation`, `material_grade`,
            `shank_length`, `shank_diameter_ratio`, `stud`) or `"MetricBolt"`
            (keys `nominal_diameter`, `thread_pitch`, `yield_stress`,
            `ultimate_tensile_stress`, `elastic_modulus`, `poissons_ratio`,
            `shank_length`, `shank_diameter_ratio`, `stud`).

    Returns:
        pyflange.bolts.MetricBolt: The bolt object.

    Raises:
        ValueError: If the bolt type is not supported.
    """
    from .bolts import StandardMetricBolt, MetricBolt

    if params['type'] == 'StandardMetricBolt':
        return StandardMetricBolt(params['designation'],
                                  params['material_grade'],
                                  shank_length = params['shank_length'],
                                  shank_diameter_ratio = params['shank_diameter_ratio'],
                                  stud = params['stud'])

    elif params['type'] == 'MetricBolt':
        return MetricBolt(nominal_diameter = params['nominal_diameter'],
                          thread_pitch = params['thread_pitch'],
                          yield_stress = params['yield_stress'],
                          ultimate_tensile_stress = params['ultimate_tensile_stress'],
                          elastic_modulus = params['elastic_modulus'],
                          poissons_ratio = params['poissons_ratio'],
                          shank_length = params['shank_length'],
                          shank_diameter_ratio = params['shank_diameter_ratio'],
                          stud = params['stud'])

    else:
        raise ValueError(f"Unsupported Bolt type '{params['type']}'")


def create_washer(params):
    """Creates a washer object from a dictionary of parameters.

    Args:
        params (dict | None): The washer parameters, or `None` for no washer.
            The `type` key selects the washer kind: `"ISOFlatWasher"` (key
            `designation`) or `"FlatWasher"` (keys `outer_diameter`,
            `inner_diameter`, `thickness`, `elastic_modulus`, `poissons_ratio`).

    Returns:
        pyflange.bolts.FlatWasher | None: The washer object, or `None`.

    Raises:
        ValueError: If the washer type is not supported.
    """
    from .bolts import ISOFlatWasher, FlatWasher

    if params is None:
        return None

    elif params['type'] == 'ISOFlatWasher':
        return ISOFlatWasher(params['designation'])

    elif params['type'] == 'FlatWasher':
        return FlatWasher(outer_diameter = params['outer_diameter'],
                          inner_diameter = params['inner_diameter'],
                          thickness = params['thickness'],
                          elastic_modulus = params['elastic_modulus'],
                          poissons_ratio = params['poissons_ratio'])

    else:
        raise ValueError(f"Unsupported Washer type '{params['type']}'")


def create_nut(**params):
    """Creates a nut object from a dictionary of parameters.

    Args:
        **params (Any): The nut parameters. The `type` key selects the nut kind:
            `"ISOHexNut"` or `"RoundNut"` (key `designation`), or `"HexNut"`
            (keys `nominal_diameter`, `thickness`, `inscribed_diameter`,
            `circumscribed_diameter`, `bearing_diameter`, `elastic_modulus`,
            `poissons_ratio`).

    Returns:
        pyflange.bolts.HexNut: The nut object.

    Raises:
        ValueError: If the nut type is not supported.
    """
    from .bolts import ISOHexNut, RoundNut, HexNut

    if params['type'] == "ISOHexNut":
        return ISOHexNut(params['designation'])

    elif params['type'] == "RoundNut":
        return RoundNut(params['designation'])

    elif params['type'] == "HexNut":
        return HexNut(nominal_diameter = params['nominal_diameter'],
                      thickness = params['thickness'],
                      inscribed_diameter = params['inscribed_diameter'],
                      circumscribed_diameter = params['circumscribed_diameter'],
                      bearing_diameter = params['bearing_diameter'],
                      elastic_modulus = params['elastic_modulus'],
                      poissons_ratio = params['poissons_ratio'])

    else:
        raise ValueError(f"Unsupported Nut type '{params['type']}'")


def create_flangesegment(**params):
    """Creates a flange segment from a dictionary of parameters.

    Given a dictionary containing the parameters for a flange segment object
    in the PyFlange JSON format (`FlangeSegment` section), returns the
    corresponding flange segment. Only L-flange segments are supported.

    The bolt preload is given either as `preload_ratio` (fraction of the bolt
    tensile yield force, `As * fy`) or as an absolute force `Fv`, but not both.

    Args:
        **params (Any): Flange segment parameters. Expected keys are 'a', 'b', 's',
            't', 'R', 'central_angle', 'Zg', 'bolt', 'preload_ratio' or 'Fv',
            'Do', 'washer', 'nut', 'gap', 'tilt_angle', 'E', 'G', 's_ratio',
            'r', 'k_shell', 'PPL', 'DFT', 'settlement_factor' and, optionally,
            'type' (default `DEFAULT_FLANGE_SEGMENT_TYPE`). Other keys are
            ignored.

    Returns:
        PolynomialLFlangeSegment: The flange segment object (possibly a
            `PolynomialLFlangeSegmentED2`).

    Raises:
        ValueError: If an unsupported type is provided, or if the preload is
            not specified exactly once.
    """
    from .flangesegments import Gap, shell_stiffness

    segment_class = _flangesegment_class(params.get('type', DEFAULT_FLANGE_SEGMENT_TYPE))
    gap = Gap(**params['gap'])
    bolt = create_bolt(**params['bolt'])

    # Bolt preload, as ratio of the yield force or as absolute value
    if ('preload_ratio' in params) == ('Fv' in params):
        raise ValueError("The flange segment preload must be specified either as 'preload_ratio' or as 'Fv'.")
    elif 'Fv' in params:
        Fv = params['Fv']
    else:
        Fv = params['preload_ratio'] * bolt.thread_cross_section.area * bolt.yield_stress

    # Shell stiffness: value, None (simplified formula) or "interp"
    k_shell = params['k_shell']
    if k_shell == "interp":
        k_shell = float(shell_stiffness(params['R'], params['s'], gap.angle))

    return segment_class(
        a = params['a'],
        b = params['b'],
        s = params['s'],
        t = params['t'],
        R = params['R'],
        central_angle = params['central_angle'],
        Zg = params['Zg'],
        bolt = bolt,
        Fv = Fv,
        Do = params['Do'],
        washer = create_washer(params.get('washer')),
        nut = create_nut(**params['nut']),
        gap = gap,
        tilt_angle = params['tilt_angle'],
        E = params['E'],
        G = params['G'],
        s_ratio = params['s_ratio'],
        r = params['r'],
        k_shell = k_shell,
        PPL = params['PPL'],
        DFT = params['DFT'],
        settlement_factor = params['settlement_factor'])


def _flangesegment_class(type_name):
    # Returns the flange segment class corresponding to a JSON type name.
    from .flangesegments import PolynomialLFlangeSegment, PolynomialLFlangeSegmentED2

    segment_classes = {
        'PolynomialLFlangeSegment': PolynomialLFlangeSegment,
        'PolynomialLFlangeSegmentED2': PolynomialLFlangeSegmentED2
    }
    if type_name not in segment_classes:
        raise ValueError(f"Unsupported FlangeSegment type: '{type_name}'")
    return segment_classes[type_name]


def create_fatiguecurve(**params):
    """Creates a fatigue curve object from a dictionary of parameters.

    Args:
        **params (Any): Fatigue curve parameters (`BoltFatigueCurve` section). The
            optional `type` key selects the curve kind:

            - `"BoltFatigueCurve"` (default): keys `BoltNominalDiameter`,
              `DS_ref`, `gamma_M`.
            - `"SingleSlopeFatigueCurve"`: keys `m`, `DS_ref`, `N_ref`.
            - `"DoubleSlopeFatigueCurve"`: keys `m1`, `m2`, `DS_ref`, `N_ref`,
              where `DS_ref` and `N_ref` define the knee point.

    Returns:
        pyflange.fatigue.FatigueCurve: The fatigue curve object.

    Raises:
        ValueError: If the curve type is not supported.
    """
    from .fatigue import BoltFatigueCurve, SingleSlopeFatigueCurve, DoubleSlopeFatigueCurve

    curve_type = params.get('type', 'BoltFatigueCurve')

    if curve_type == 'BoltFatigueCurve':
        return BoltFatigueCurve(params['BoltNominalDiameter'],
                                DS_ref = params['DS_ref'],
                                gamma_M = params['gamma_M'])

    elif curve_type == 'SingleSlopeFatigueCurve':
        return SingleSlopeFatigueCurve(params['m'], params['DS_ref'], params['N_ref'])

    elif curve_type == 'DoubleSlopeFatigueCurve':
        return DoubleSlopeFatigueCurve(params['m1'], params['m2'], params['DS_ref'], params['N_ref'])

    else:
        raise ValueError(f"Unsupported fatigue curve type '{curve_type}'")


def create_markovmatrix(**params):
    """Creates a Markov matrix object from a dictionary of parameters.

    Args:
        **params (Any): Markov matrix parameters (`MarkovMatrix` section).
            Expected keys are 'range', 'mean', 'cycles', 'duration'.

    Returns:
        pyflange.fatigue.MarkovMatrix: The Markov matrix object.
    """
    from .fatigue import MarkovMatrix
    return MarkovMatrix(range    = np.array(params['range'], dtype=float),
                        mean     = np.array(params['mean'], dtype=float),
                        cycles   = np.array(params['cycles'], dtype=float),
                        duration = params['duration'])


def create_fatiguecase(data, SMF=None):
    """Creates a BoltFatigueAnalysis object from a case dictionary.

    Args:
        data (dict): A case dictionary, in any supported format version.
        SMF (float, optional): Stress multiplication factor. If omitted, the
            `SMF` value of the case is used, or 1.0 if the case has none.

    Returns:
        pyflange.fatigue.BoltFatigueAnalysis: The fatigue analysis object.
    """
    from .fatigue import BoltFatigueAnalysis

    data = upgrade_case(data)
    return BoltFatigueAnalysis(create_flangesegment(**data['FlangeSegment']),
                               create_markovmatrix(**data['MarkovMatrix']),
                               create_fatiguecurve(**data['BoltFatigueCurve']),
                               allowable_damage = data.get('AllowableDamage', data.get('allowable_damage', 1.0)),
                               SMF = SMF if SMF is not None else data.get('SMF', 1.0))


def load_fatiguecase(file_path, SMF=None):
    """Loads a BoltFatigueAnalysis object from a JSON file.

    Given the path of a JSON file containing a flange fatigue case in the
    PyFlange JSON format, returns the corresponding BoltFatigueAnalysis object.

    Args:
        file_path (str): The path to the JSON file.
        SMF (float, optional): Stress multiplication factor. If omitted, the
            `SMF` value of the case is used, or 1.0 if the case has none.

    Returns:
        pyflange.fatigue.BoltFatigueAnalysis: The loaded fatigue analysis object.
    """
    return create_fatiguecase(load_inputfile(file_path), SMF=SMF)


def create_distribution(**params):
    """Creates a scipy distribution from a dictionary of parameters.

    Args:
        **params (Any): Distribution parameters. The `Distribution` key selects
            the distribution kind:

            - `"Normal"`: keys `Mean`, `CoV`.
            - `"LogNormal"`: keys `Mean`, `CoV`.
            - `"Weibull"`: keys `Shape`, `Scale`.
            - `"Beta"`: keys `Alpha`, `Beta` and optionally `x0`, `x1`
              (bounds of the definition interval, default 0 and 1).

    Returns:
        scipy.stats.rv_frozen: A frozen scipy continuous distribution.

    Raises:
        ValueError: If an unsupported distribution type is provided.
    """
    import scipy.stats
    from math import log, exp, sqrt

    if params['Distribution'] == "Normal":
        return scipy.stats.norm(params['Mean'], params['CoV'] * params['Mean'])

    elif params['Distribution'] == "LogNormal":
        shape = sqrt( log(params['CoV']**2 + 1) )
        scale = exp(log(params['Mean']) - shape**2 / 2)
        return scipy.stats.lognorm(s=shape, loc=0, scale=scale)

    elif params['Distribution'] == "Weibull":
        return scipy.stats.weibull_min(params['Shape'], scale=params['Scale'])

    elif params['Distribution'] == "Beta":
        x0 = params.get('x0', 0)
        x1 = params.get('x1', 1)
        return scipy.stats.beta(params['Alpha'], params['Beta'], loc=x0, scale=x1-x0)

    else:
        raise ValueError(f"Unsupported probability distribution '{params['Distribution']}'")


def create_sampler(**params):
    """Creates a pyflange sampler from a dictionary of parameters.

    Args:
        **params (Any): Sampler parameters. All the distributions supported by
            `create_distribution` are accepted, plus `"Product"`, which
            samples the product of two independent variables described by
            the `Dist1` and `Dist2` keys.

    Returns:
        Generator: A sampler yielding random values.

    Raises:
        ValueError: If an unsupported distribution type is provided.
    """
    from . import stats

    if params['Distribution'] == "Product":
        return stats.product_sampler(create_distribution(**params['Dist1']),
                                     create_distribution(**params['Dist2']))
    else:
        return stats.sampler(create_distribution(**params))


def create_gap_sampler(flange_diameter, **params):
    """Creates a random gap sampler from a dictionary of parameters.

    Args:
        flange_diameter (float): The outer diameter of the flange.
        **params (Any): Gap sampler parameters (`Stats.Gap` section). The
            `JointDistribution` key selects the kind of sampler:

            - `"IEC 61400-6:2020 AMD1"`: keys `FlatnessTolerance`, `Angle`
              (distribution), `ShapeFactor` (distribution).
            - `"IEC 61400-6:2020 ED2"`: key `FlatnessTolerance`.
            - `"ClaytonCopula"`: keys `Theta`, `Height` (distribution),
              `Angle` (distribution), `ShapeFactor` (distribution).

    Returns:
        Generator: A sampler yielding random `pyflange.flangesegments.Gap` objects.

    Raises:
        ValueError: If an unsupported joint distribution is provided.
    """
    from . import stats

    if params['JointDistribution'] == "IEC 61400-6:2020 AMD1":
        return stats.standard_gap_sampler(
                flange_diameter = flange_diameter,
                flange_flatness_tolerance = params['FlatnessTolerance'],
                gap_angle_sampler = create_sampler(**params['Angle']),
                gap_shape_factor_sampler = create_sampler(**params['ShapeFactor']))

    elif params['JointDistribution'] == "IEC 61400-6:2020 ED2":
        return stats.standard_gap_sampler_ED2(
                flange_diameter = flange_diameter,
                flange_flatness_tolerance = params['FlatnessTolerance'])

    elif params['JointDistribution'] == "ClaytonCopula":
        return stats.clayton_gap_sampler(
                theta = params['Theta'],
                height_marginal_dist = create_distribution(**params['Height']),
                angle_marginal_dist  = create_distribution(**params['Angle']),
                shape_factor_sampler = create_sampler(**params['ShapeFactor']))

    else:
        raise ValueError(f"Unsupported gap joint distribution '{params['JointDistribution']}'")


def create_fatiguecase_sampler(fseg, markov_matrix, **params):
    """Creates a sampler of random fatigue cases.

    The random flange segments are of the same type as `fseg` and share its
    deterministic parameters. Their shell stiffness is always interpolated
    from the FEA data (see `pyflange.flangesegments.shell_stiffness`) as a
    function of the random gap angle.

    Args:
        fseg (PolynomialLFlangeSegment): The base flange segment object.
        markov_matrix (MarkovMatrix): The base Markov matrix object.
        **params (Any): Statistical parameters (`Stats` section, current format
            version), for preload, plastic loss, DFT, gap, tilt, load range
            factor, fatigue class factor, allowable damage and stress
            multiplication factor.

    Returns:
        Generator: A generator that yields randomized BoltFatigueAnalysis objects.

    Raises:
        ValueError: If an unsupported distribution type is provided for any
            parameter.
        KeyError: If a required parameter is missing.
    """
    from . import stats

    # Preload sampler, with mean value given as ratio of the yield force or as absolute value
    preload = params['Preload']
    if 'MeanRatio' in preload:
        mean_preload = preload['MeanRatio'] * fseg.bolt.thread_cross_section.area * fseg.bolt.yield_stress
    else:
        mean_preload = preload['Mean']

    if preload['Distribution'] == "Normal":
        preload_samp = stats.norm_sampler(mean_preload, preload['CoV'])
    elif preload['Distribution'] == "LogNormal":
        preload_samp = stats.lognorm_sampler(mean_preload, preload['CoV'])
    else:
        raise ValueError(f"Unsupported preload distribution '{preload['Distribution']}'")

    # Flange-segment sampler
    fseg_samp = stats.standard_PolynomialLFlangeSegment_sampler(
        a = fseg.a, b = fseg.b, s = fseg.s, t = fseg.t, R = fseg.R,
        central_angle = fseg.central_angle,
        Zg = fseg.Zg,
        bolt = fseg.bolt, Do = fseg.Do, washer = fseg.washer, nut = fseg.nut,
        preload_sampler = preload_samp,
        ppl_sampler = create_sampler(**params['PlasticLoss']),
        dft_sampler = create_sampler(**params['DFT']),
        gap_sampler = create_gap_sampler(2*fseg.R, **params['Gap']),
        tilt_sampler = create_sampler(**params['Tilt']),
        E = fseg.E, G = fseg.G, s_ratio = fseg.s_ratio, r = fseg.r,
        k_shell = "interp",
        settlement_factor = fseg.settlement_factor,
        segment_class = type(fseg))

    # Uncertainty in load ranges
    load_range = params['LoadRangeFactor']
    if load_range['Distribution'] != "LogNormal":
        raise ValueError(f"Unsupported load range factor distribution '{load_range['Distribution']}'")
    markov_matrix_samp = stats.standard_markov_matrix_sampler(markov_matrix, load_range['Mean'], load_range['CoV'])

    # Uncertainty in fatigue capacity
    fatigue_class = params['FatigueClassFactor']
    if fatigue_class['Distribution'] != "Normal":
        raise ValueError(f"Unsupported fatigue class factor distribution '{fatigue_class['Distribution']}'")
    if 'DS_ref_mean' not in fatigue_class:
        raise KeyError("Missing 'DS_ref_mean' in the 'FatigueClassFactor' parameters.")
    fatigue_curve_samp = stats.standard_bolt_fatigue_curve_sampler(
            fseg.bolt.nominal_diameter, fatigue_class['Mean'], fatigue_class['CoV'], fatigue_class['DS_ref_mean'])

    # Stress multiplication factor: sampler or constant value
    smf = params['StressMultiplicationFactor']
    smf = create_sampler(**smf) if isinstance(smf, dict) else float(smf)

    return stats.fatigue_case_sampler(fseg_samp,
                                      markov_matrix_samp,
                                      fatigue_curve_samp,
                                      create_sampler(**params['AllowableDamage']),
                                      SMF = smf)


def load_fatiguecase_sampler(file_path):
    """Loads a fatigue case sampler from a JSON file.

    Args:
        file_path (str): The path to the JSON file containing the fatigue case
            and its `Stats` section.

    Returns:
        Generator: A generator that yields randomized BoltFatigueAnalysis objects.
    """
    data = upgrade_case(load_inputfile(file_path))
    fseg = create_flangesegment(**data['FlangeSegment'])
    markov_matrix = create_markovmatrix(**data['MarkovMatrix'])
    return create_fatiguecase_sampler(fseg, markov_matrix, **data['Stats'])





# -----------------------------------------------------------------------------
#   MONTE CARLO SIMULATION
# -----------------------------------------------------------------------------

def sample_fatiguecase(data, size, seed=None, progress=None):
    """Runs a Monte Carlo simulation of a fatigue case.

    Args:
        data (dict): A case dictionary, including the `Stats` section.
        size (int): Number of realizations to generate.
        seed (int, optional): Seed of the random number generator. If given,
            the numpy global random generator is re-seeded, so that the
            simulation is reproducible. If omitted, the current random state
            is used.
        progress (callable, optional): A function called as `progress(i)`
            after the generation of the i-th realization (1-based).

    Returns:
        pandas.DataFrame: A table with one row per realization and the
            columns listed in `SAMPLE_COLUMNS`. Lengths are in m, forces in N,
            angles in rad, stresses in Pa, fatigue life in the time units of
            the Markov matrix duration. `FatigueClass` is the knee-point
            stress range of the sampled bolt S-N curve and
            `MomentRangeCoefficient` the factor applied to the load ranges.
    """
    import pandas as pd

    data = upgrade_case(data)
    if seed is not None:
        np.random.seed(seed)

    fseg = create_flangesegment(**data['FlangeSegment'])
    markov_matrix = create_markovmatrix(**data['MarkovMatrix'])
    fcase_samp = create_fatiguecase_sampler(fseg, markov_matrix, **data['Stats'])

    # Reference load range, to evaluate the sampled load range coefficient
    nonzero = np.flatnonzero(markov_matrix.range)
    i_ref = nonzero[0] if len(nonzero) else 0
    ref_range = markov_matrix.range[i_ref]

    records = []
    for i in range(int(size)):
        rnd = next(fcase_samp)
        records.append({
            'BoltPreload'           : rnd.fseg.bolt_preload,
            'GapShapeFactor'        : rnd.fseg.gap.shape_factor,
            'TiltAngle'             : rnd.fseg.tilt_angle,
            'GapAngle'              : rnd.fseg.gap.angle,
            'GapHeight'             : rnd.fseg.gap.height,
            'FatigueClass'          : rnd.fatigue_curve.curves[0].DS_ref,
            'MomentRangeCoefficient': rnd.flange_mkvm.range[i_ref] / ref_range,
            'Damage'                : rnd.damage,
            'AllowableDamage'       : rnd.allowable_damage,
            'FatigueLife'           : rnd.fatigue_life,
            'SMF'                   : rnd.SMF,
            'Utilization'           : rnd.damage / rnd.allowable_damage
        })
        if progress:
            progress(i+1)

    return pd.DataFrame.from_records(records, columns=SAMPLE_COLUMNS).astype(float)


def create_realization(data, record):
    """Rebuilds a single Monte Carlo realization.

    Given the case that generated a sample and one row of the table returned
    by `sample_fatiguecase` (or of the corresponding CSV file), rebuilds the
    random fatigue case that produced that row.

    Args:
        data (dict): The case dictionary that generated the sample.
        record (Mapping): One realization, with the keys listed in
            `SAMPLE_COLUMNS` (e.g. a `pandas.Series` row).

    Returns:
        pyflange.fatigue.BoltFatigueAnalysis: The rebuilt fatigue case.
    """
    from .fatigue import BoltFatigueAnalysis, DoubleSlopeFatigueCurve
    from .stats import _shell_stiffness_interpolator

    data = upgrade_case(data)
    fseg_params = data['FlangeSegment']

    # The sampled preload already includes all the losses
    fseg_params.pop('preload_ratio', None)
    fseg_params['Fv'] = float(record['BoltPreload'])
    fseg_params['PPL'] = 0.0
    fseg_params['DFT'] = 0.0

    # Sampled gap, tilt and corresponding shell stiffness
    gap_angle = float(record['GapAngle'])
    fseg_params['gap'] = {
        'height': float(record['GapHeight']),
        'angle': gap_angle,
        'shape_factor': float(record['GapShapeFactor'])
    }
    fseg_params['tilt_angle'] = float(record['TiltAngle'])
    fseg_params['k_shell'] = float(_shell_stiffness_interpolator(fseg_params['R'], fseg_params['s'])(gap_angle))

    # Sampled loads and bolt S-N curve (FatigueClass is the knee-point stress range)
    mkvm_params = data['MarkovMatrix']
    mkvm_params['range'] = list(np.array(mkvm_params['range'], dtype=float) * float(record['MomentRangeCoefficient']))
    fatigue_curve = DoubleSlopeFatigueCurve(3, 5, float(record['FatigueClass']), 2.0e6)

    return BoltFatigueAnalysis(create_flangesegment(**fseg_params),
                               create_markovmatrix(**mkvm_params),
                               fatigue_curve,
                               allowable_damage = float(record['AllowableDamage']),
                               SMF = float(record['SMF']))





# -----------------------------------------------------------------------------
#   EXPORT DATA
# -----------------------------------------------------------------------------

def bolt_to_dict(bolt):
    """Converts a bolt object to a dictionary in the JSON format.

    Args:
        bolt (pyflange.bolts.MetricBolt): The bolt object.

    Returns:
        dict: The bolt parameters, with type `"MetricBolt"`.

    Raises:
        ValueError: If the bolt type is not supported.
    """
    from .bolts import MetricBolt

    if not isinstance(bolt, MetricBolt):
        raise ValueError(f"Unknown bolt type: '{type(bolt)}'")

    return {
        "type": "MetricBolt",
        "nominal_diameter": float(bolt.nominal_diameter),
        "thread_pitch": float(bolt.thread_pitch),
        "yield_stress": float(bolt.yield_stress),
        "ultimate_tensile_stress": float(bolt.ultimate_tensile_stress),
        "elastic_modulus": float(bolt.elastic_modulus),
        "poissons_ratio": float(bolt.poissons_ratio),
        "shank_length": float(bolt.shank_length),
        "shank_diameter_ratio": float(bolt.shank_diameter_ratio),
        "stud": bool(bolt.stud)
    }


def washer_to_dict(washer):
    """Converts a washer object to a dictionary in the JSON format.

    Args:
        washer (pyflange.bolts.FlatWasher | None): The washer object.

    Returns:
        dict | None: The washer parameters, with type `"FlatWasher"`, or
            `None` if `washer` is `None`.

    Raises:
        ValueError: If the washer type is not supported.
    """
    from .bolts import FlatWasher

    if washer is None:
        return None

    if not isinstance(washer, FlatWasher):
        raise ValueError(f"Unknown washer type: '{type(washer)}'")

    return {
        "type": "FlatWasher",
        "outer_diameter": float(washer.outer_diameter),
        "inner_diameter": float(washer.inner_diameter),
        "thickness": float(washer.thickness),
        "elastic_modulus": float(washer.elastic_modulus),
        "poissons_ratio": float(washer.poissons_ratio)
    }


def nut_to_dict(nut):
    """Converts a nut object to a dictionary in the JSON format.

    Args:
        nut (pyflange.bolts.HexNut): The nut object.

    Returns:
        dict: The nut parameters, with type `"HexNut"`.

    Raises:
        ValueError: If the nut type is not supported.
    """
    from .bolts import HexNut

    if not isinstance(nut, HexNut):
        raise ValueError(f"Unknown nut type: '{type(nut)}'")

    return {
        "type": "HexNut",
        "nominal_diameter": float(nut.nominal_diameter),
        "thickness": float(nut.thickness),
        "inscribed_diameter": float(nut.inscribed_diameter),
        "circumscribed_diameter": float(nut.circumscribed_diameter),
        "bearing_diameter": float(nut.bearing_diameter),
        "elastic_modulus": float(nut.elastic_modulus),
        "poissons_ratio": float(nut.poissons_ratio)
    }


def gap_to_dict(gap):
    """Converts a gap object to a dictionary in the JSON format.

    Args:
        gap (pyflange.flangesegments.Gap): The gap object.

    Returns:
        dict: The gap parameters.

    Raises:
        ValueError: If `gap` is not a Gap object.
    """
    from .flangesegments import Gap

    if not isinstance(gap, Gap):
        raise ValueError("Invalid gap object.")

    return {
        "height": float(gap.height),
        "angle": float(gap.angle),
        "shape_factor": float(gap.shape_factor)
    }


def flangesegment_to_dict(fseg, verbose=False):
    """Converts a flange segment to a dictionary in the JSON format.

    Args:
        fseg (PolynomialLFlangeSegment): The flange segment object.
        verbose (bool, optional): If True, the dictionary also contains
            a report of derived and intermediate quantities (preload after
            losses, polynomial reference points and coefficients, transfer
            functions, stiffnesses, etc.). These extra keys are ignored when
            the dictionary is loaded back. Defaults to False.

    Returns:
        dict: The flange segment parameters.

    Raises:
        ValueError: If the flange segment type is not supported.
    """
    from .utils import read_data_log

    type_name = type(fseg).__name__
    if type_name not in FLANGE_SEGMENT_TYPES:
        raise ValueError(f"Unknown flange segment type: '{type(fseg)}'")

    k_shell = fseg.k_shell
    d = {
        "type": type_name,

        "a": fseg.a,        # inner face to bolt center
        "b": fseg.b,        # bolt center to shell centerline
        "s": fseg.s,
        "t": fseg.t,
        "R": fseg.R,
        "central_angle": fseg.central_angle,

        "Zg": fseg.Zg,

        "bolt"  : bolt_to_dict(fseg.bolt),
        "preload_ratio" : fseg.Fv / (fseg.bolt.thread_cross_section.area * fseg.bolt.yield_stress),
        "Do"    : fseg.Do,
        "washer": washer_to_dict(fseg.washer),
        "nut"   : nut_to_dict(fseg.nut),

        "gap"       : gap_to_dict(fseg.gap),
        "tilt_angle": fseg.tilt_angle,

        "E": fseg.E,
        "G": fseg.G,

        "s_ratio": fseg.s_ratio,
        "r"      : fseg.r,
        "k_shell": None if k_shell is None else float(k_shell),

        "PPL": fseg.PPL,
        "DFT": fseg.DFT,
        "settlement_factor": fseg.settlement_factor
    }

    if verbose:
        d["BoltPreload"] = fseg.bolt_preload
        d["Point1"] = {
            "Z" : fseg.shell_force_at_rest,
            "Fs": fseg.bolt_force_at_rest,
            "Ms": fseg.bolt_moment_at_rest
        }
        d["Point2"] = {
            "Z" : fseg.shell_force_at_tensile_ULS,
            "Fs": fseg.bolt_force_at_tensile_ULS,
            "Ms": fseg.bolt_moment_at_tensile_ULS
        }
        d["Point3"] = {
            "Z" : fseg.shell_force_at_small_displacement,
            "Fs": fseg.bolt_force_at_small_displacement,
            "Ms": fseg.bolt_moment_at_small_displacement
        }
        Z4 = fseg._compressive_force_polynomial.domain[0]
        d["Point4"] = {
            "Z" : Z4,
            "Fs": fseg.bolt_axial_force(Z4),
            "Ms": fseg.bolt_bending_moment(Z4)
        }
        d["AxialForcePolinomialCoeff"] = {
            "TensileZone": list(fseg._tensile_force_polynomial.coef),
            "CompressiveZone": list(fseg._compressive_force_polynomial.coef),
            "ClosedGapZone": [d["Point4"]["Fs"], 0, 0]
        }
        d["BendingMomentPolinomialCoeff"] = {
            "TensileZone": list(fseg._tensile_moment_polynomial.coef),
            "CompressiveZone": list(fseg._compressive_moment_polynomial.coef),
            "ClosedGapZone": [d["Point4"]["Ms"], 0, 0]
        }
        Z = np.linspace(1.05*fseg._compressive_force_polynomial.domain[0], 1.05*fseg.shell_force_at_tensile_ULS, 100)
        d["TransferFunctions"] = {
            "Z" : list(Z),
            "Fs": list(fseg.bolt_axial_force(Z)),
            "Ms": list(fseg.bolt_bending_moment(Z))
        }
        d["bolt_As"] = fseg.bolt.thread_cross_section.area
        d["bolt_axial_stiffness"] = fseg._bolt_axial_stiffness
        d["a_prime"] = fseg.b / (fseg._prying_lever_ratio - 1)
        d["Z0"] = fseg._ideal_shell_force_at_tensile_ULS
        d["bolt.bending_stiffness"] = fseg._bolt_bending_stiffness
        d["clamped_parts_stiffness"] = fseg._flange_axial_stiffness
        d["gap.stiffness"] = fseg._gap_stiffness
        d["DZ_gap_inclination"] = fseg._tilt_neutralization_shell_force
        d["DZ_gap"] = fseg._parallel_gap_neutralization_shell_force
        d["DZ_gap_tot"] = fseg._total_gap_neutralization_shell_force
        d["stiffness_correction_factor"] = fseg._stiffness_correction_factor
        d["polynomial_initial_slope"] = fseg._polynomial_initial_slope
        d["true_force_initial_slope"] = fseg._tensile_force_polynomial.deriv()(d["Point1"]["Z"])
        d["true_moment_initial_slope"] = fseg._tensile_moment_polynomial.deriv()(d["Point1"]["Z"])
        d["gap_k_mean"] = read_data_log(fseg, "k_mean")
        d["gap_COV"] = read_data_log(fseg, "COV_k")
        d["gap_k_fac"] = read_data_log(fseg, "k_fac")
        d["gap_k_shell"] = fseg._k_shell
        d["gap_stiffness_correction_factor"] = read_data_log(fseg, "f_tot")
        d["flange_A_cf"] = read_data_log(fseg, "A_cf")
        d["flange_I_cf"] = read_data_log(fseg, "I_cf")
        d["gap_k_fl"] = read_data_log(fseg, "k_fl")
        d["u"] = read_data_log(fseg, "u")
        d["k_seg"] = read_data_log(fseg, "k_seg")
        d["a_star"] = fseg._a_star
        d["I_tg"] = read_data_log(fseg, "I_tg")

    return _plain(d)


def fatiguecurve_to_dict(sncurve):
    """Converts a fatigue curve object to a dictionary in the JSON format.

    Args:
        sncurve (pyflange.fatigue.FatigueCurve): A `BoltFatigueCurve`,
            `SingleSlopeFatigueCurve` or `DoubleSlopeFatigueCurve` object.

    Returns:
        dict: The fatigue curve parameters. For a `BoltFatigueCurve`, the
            knee point of the design curve is reported too, as `m1`, `m2`,
            `DS12` and `N_ref`.

    Raises:
        ValueError: If the fatigue curve type is not supported.
    """
    from .fatigue import BoltFatigueCurve, SingleSlopeFatigueCurve, DoubleSlopeFatigueCurve

    if isinstance(sncurve, BoltFatigueCurve):
        return _plain({
            "type": "BoltFatigueCurve",
            "BoltNominalDiameter": sncurve.diameter,
            "DS_ref": sncurve.DS_ref,
            "gamma_M": sncurve.gamma_M,
            "m1": sncurve.curves[0].m,
            "m2": sncurve.curves[1].m,
            "DS12": sncurve.curves[0].DS_ref,
            "N_ref": sncurve.curves[0].N_ref
        })

    elif isinstance(sncurve, SingleSlopeFatigueCurve):
        return _plain({
            "type": "SingleSlopeFatigueCurve",
            "m": sncurve.m,
            "DS_ref": sncurve.DS_ref,
            "N_ref": sncurve.N_ref
        })

    elif isinstance(sncurve, DoubleSlopeFatigueCurve):
        return _plain({
            "type": "DoubleSlopeFatigueCurve",
            "m1": sncurve.curves[0].m,
            "m2": sncurve.curves[1].m,
            "DS_ref": sncurve.curves[0].DS_ref,
            "N_ref": sncurve.curves[0].N_ref
        })

    else:
        raise ValueError(f"Unknown fatigue curve object '{type(sncurve)}'")


def markovmatrix_to_dict(mkvm):
    """Converts a Markov matrix object to a dictionary in the JSON format.

    Args:
        mkvm (pyflange.fatigue.MarkovMatrix): The Markov matrix object.

    Returns:
        dict: The Markov matrix parameters.

    Raises:
        ValueError: If `mkvm` is not a MarkovMatrix object.
    """
    from .fatigue import MarkovMatrix

    if not isinstance(mkvm, MarkovMatrix):
        raise ValueError("Invalid Markov Matrx object.")

    return _plain({
        "range"    : list(mkvm.range),
        "mean"     : list(mkvm.mean),
        "cycles"   : list(mkvm.cycles),
        "duration" : mkvm.duration
    })


def fatiguecase_to_dict(fcase, verbose=False):
    """Converts a fatigue analysis to a dictionary in the JSON format.

    The returned dictionary contains the case input (it can be loaded back
    with `create_fatiguecase`) and the analysis results.

    Args:
        fcase (pyflange.fatigue.BoltFatigueAnalysis): The fatigue analysis.
        verbose (bool, optional): If True, the flange segment report (see
            `flangesegment_to_dict`), the bolt stress Markov matrix and the
            fatigue life are included too. Defaults to False.

    Returns:
        dict: The fatigue case data and results (`damage`, `allowable_damage`,
            `SMF` and, if verbose, `FatigueLife`).

    Raises:
        ValueError: If `fcase` is not a BoltFatigueAnalysis object.
    """
    from .fatigue import BoltFatigueAnalysis

    if not isinstance(fcase, BoltFatigueAnalysis):
        raise ValueError("Invalid fatigue case object")

    d = {
        "format_version": FORMAT_VERSION,
        "FlangeSegment": flangesegment_to_dict(fcase.fseg, verbose),
        "BoltFatigueCurve": fatiguecurve_to_dict(fcase.fatigue_curve),
        "MarkovMatrix": markovmatrix_to_dict(fcase.flange_mkvm),
        "damage": fcase.damage,
        "allowable_damage": fcase.allowable_damage,
        "SMF": fcase.SMF
    }

    if verbose:
        d["BoltStressMarkovMatrix"] = markovmatrix_to_dict(fcase.bolt_mkvm)
        d["FatigueLife"] = fcase.fatigue_life

    return _plain(d)


def dict_to_json(d, file_path):
    """Writes a dictionary to a JSON file.

    Numpy scalars and arrays are converted to plain Python numbers and lists.

    Args:
        d (dict): The dictionary to be written.
        file_path (str): The path to the JSON file to be created.
    """
    with open(file_path, "w") as file:
        json.dump(_plain(d), file, indent=4)


def dict_to_excel(d, file_path):
    """Writes a dictionary to an Excel file.

    This utility function takes a dictionary and an Excel file name as inputs,
    flattens the dictionary into a list of chained-key-value pairs, and
    writes it to an Excel file in a two-column table, where the first column
    contains the keys and the second column contains the corresponding values.
    Each value cell is also given an Excel name equal to its key, so that it
    can be referenced by name from other workbooks.

    Example:
        Given the following dictionary:
            d = {
                "a": 10,
                "b": {
                    "c": 20,
                    "d": 30
                }
            }

        This function will write it to Excel in the following format:
            | KEY | VALUE |
            |-----|-------|
            | a   | 10    |
            | b.c | 20    |
            | b.d | 30    |

    Args:
        d (dict): The dictionary to be written.
        file_path (str): The path to the Excel file to be created.
    """
    import openpyxl
    from openpyxl.workbook.defined_name import DefinedName

    def iterate_flat_dict(d):
        # Recursively flattens a nested dictionary into (chained key, value) pairs.
        for key, value in d.items():
            if isinstance(value, dict):
                for sub_key, sub_value in iterate_flat_dict(value):
                    yield key + "." + sub_key, sub_value
            elif isinstance(value, list):
                yield key, ", ".join([str(item) for item in value])
            else:
                yield key, value

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "CaseData"

    ws.append(["KEY", "VALUE"])
    i = 2
    for key, value in iterate_flat_dict(_plain(d)):
        ws.append([key, value])
        ws.defined_names.add(DefinedName(key, attr_text=f"{ws.title}!$B${i}"))
        i += 1

    wb.save(file_path)


def _plain(value):
    # Recursively converts numpy scalars and arrays to plain python types.
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    elif isinstance(value, np.ndarray) and value.ndim == 0:
        return _plain(value.item())
    elif isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    elif isinstance(value, np.bool_):
        return bool(value)
    elif isinstance(value, np.integer):
        return int(value)
    elif isinstance(value, np.floating):
        return float(value)
    else:
        return value
