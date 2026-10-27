PyFlange JSON case format
=========================

A *fatigue case* describes a single-bolt flange segment, its fatigue loads,
its bolt S-N curve and, optionally, the statistical parameters needed for a
Monte Carlo simulation. PyFlange reads and writes fatigue cases as JSON files
in the format described on this page.

The format is used by the [command-line interface](cli.md) and by the
[`pyflange.io`](apis/io.md) module:

```python
from pyflange.io import load_fatiguecase, load_fatiguecase_sampler

fcase = load_fatiguecase("case.json")           # deterministic analysis
print(fcase.damage, fcase.fatigue_life)

samp = load_fatiguecase_sampler("case.json")    # random fatigue cases
rnd_fcase = next(samp)
```

A complete, working example is available here:
[fatigue-case.json](examples/fatigue-case.json).

> **Scope**: the format currently supports **L-flange segments only**
> (`PolynomialLFlangeSegment` and `PolynomialLFlangeSegmentED2`).
> T-flange segments (`PolynomialTFlangeSegment`) are not supported yet.


Format versions
---------------

The format version is stored in the top-level `format_version` field. The
current version is **0.2.0**.

| Version | Notes |
|---------|-------|
| `0.1.0` | Original format of the Forrester tool. Files of this version have **no** `format_version` field. |
| `0.2.0` | Current version. See [changes from 0.1.0](#changes-from-version-010). |

Files of any supported version can be read: files without `format_version`
are assumed to be version 0.1.0 and are converted to the current version on
load (see `pyflange.io.upgrade_case`), so that they keep producing the same
results. PyFlange always writes the current version.


Conventions
-----------

- **Units**: all values are in SI units: lengths in m, forces in N, moments
  in N·m, stresses and elastic moduli in Pa, angles in **radians**. The only
  exception is the Markov matrix `duration`, whose unit is free and
  determines the unit of the resulting fatigue life (typically years).
- **No hidden defaults**: all the parameters marked as *required* must be
  given explicitly. In particular, PyFlange does not provide default values
  for any statistical parameter.
- **Unknown keys are ignored**. This allows adding notes or extra data to a
  case file, and allows reading back the result files written by PyFlange
  (see [Result files](#result-files)), which contain extra report keys.
- Keys are case-sensitive.


Top-level structure
-------------------

```json
{
    "format_version": "0.2.0",
    "CaseID": "My case",
    "FlangeSegment": { ... },
    "BoltFatigueCurve": { ... },
    "MarkovMatrix": { ... },
    "AllowableDamage": 1.0,
    "SMF": 1.0,
    "Stats": { ... }
}
```

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `format_version` | string | yes (0.2.0) | Format version. If missing, the file is read as version 0.1.0. |
| `CaseID` | any | no | Case identifier. Not used by PyFlange. |
| `FlangeSegment` | object | yes | The flange segment. See [FlangeSegment](#flangesegment). |
| `BoltFatigueCurve` | object | yes | The bolt S-N curve of the deterministic analysis. See [BoltFatigueCurve](#boltfatiguecurve). |
| `MarkovMatrix` | object | yes | The flange fatigue loads. See [MarkovMatrix](#markovmatrix). |
| `AllowableDamage` | number | no | Allowable damage of the deterministic analysis. If omitted, 1.0. |
| `SMF` | number | no | Stress multiplication factor of the deterministic analysis, applied to the bolt stress ranges. If omitted, 1.0. The command-line option `--smf` overrides it. |
| `Stats` | object | for Monte Carlo only | Statistical parameters. See [Stats](#stats). |


FlangeSegment
-------------

```json
"FlangeSegment": {
    "type": "PolynomialLFlangeSegment",
    "a": 0.150, "b": 0.122, "s": 0.054, "t": 0.172, "R": 4.000,
    "central_angle": 0.0402768,
    "Zg": -115667,
    "bolt": { "type": "StandardMetricBolt", "designation": "M80", "material_grade": "10.9",
              "shank_length": 0.160, "shank_diameter_ratio": "pitch", "stud": true },
    "preload_ratio": 0.75301,
    "Do": 0.086,
    "washer": null,
    "nut": { "type": "RoundNut", "designation": "M80" },
    "gap": { "height": 0.00089, "angle": 0.5236, "shape_factor": 1.0 },
    "tilt_angle": 0.001745,
    "E": 210e9, "G": 80.77e9,
    "s_ratio": 1.0, "r": 0.01,
    "k_shell": "interp",
    "PPL": 0.05, "DFT": 300e-6,
    "settlement_factor": 0.5
}
```

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `type` | string | - | `"PolynomialLFlangeSegment"` (IEC 61400-6 AMD1 model) or `"PolynomialLFlangeSegmentED2"` (IEC 61400-6 ED2 model). Optional, default `"PolynomialLFlangeSegment"`. |
| `a` | number | m | Distance between the inner face of the flange and the center of the bolt hole. |
| `b` | number | m | Distance between the center of the bolt hole and the center-line of the shell. |
| `s` | number | m | Shell thickness. |
| `t` | number | m | Flange thickness. |
| `R` | number | m | Shell outer radius. |
| `central_angle` | number | rad | Angle subtended by the flange segment (2π / number of bolts). |
| `Zg` | number | N | Shell force per segment at rest (dead weight divided by the number of bolts). Negative if compressive. |
| `bolt` | object | - | The bolt. See [Bolts](#bolts). |
| `preload_ratio` | number | - | Bolt preload before losses, as a fraction of the bolt tensile yield force `As·fy`, where `As` is the area of the bolt thread cross-section and `fy` the bolt yield stress. Exactly one of `preload_ratio` and `Fv` is required. |
| `Fv` | number | N | Bolt preload before losses, as absolute value. Exactly one of `preload_ratio` and `Fv` is required. |
| `Do` | number | m | Bolt hole diameter. Standard values are given by `pyflange.bolts.standard_hole_diameter`. |
| `washer` | object or null | - | The washer, or `null` (or omitted) for no washer. See [Washers](#washers). |
| `nut` | object | - | The nut. See [Nuts](#nuts). |
| `gap` | object | - | The flange gap. See [Gap](#gap). |
| `tilt_angle` | number | rad | Flange tilt angle. |
| `E` | number | Pa | Young's modulus of the flange. |
| `G` | number | Pa | Shear modulus of the flange. |
| `s_ratio` | number | - | Ratio of the bottom shell thickness over `s`. |
| `r` | number | m | Rounding radius between flange and shell. |
| `k_shell` | number, null or `"interp"` | N/m/m | Shell stiffness. A number is used as it is; `null` selects the simplified formula of the model; `"interp"` interpolates the value from FEA data (see `pyflange.flangesegments.shell_stiffness`) for the given `R`, `s` and gap angle. |
| `PPL` | number | - | Preload loss due to plastic deformations, as a fraction of `Fv` (e.g. 0.05 for 5%). |
| `DFT` | number | m | Dry film thickness of the coatings, used to evaluate the settlement preload losses. |
| `settlement_factor` | number | - | Factor reducing the settlement preload losses (e.g. 1.0 for torque tightening, 0.5 for tension tightening). |

All keys are required, except `type`, `washer` and one of `preload_ratio` / `Fv`.

### Bolts

The `type` key selects one of the following bolt kinds.

**`StandardMetricBolt`**: a bolt with standard ISO dimensions and material
(see `pyflange.bolts.StandardMetricBolt`).

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `designation` | string | - | Metric thread designation, e.g. `"M80"`. |
| `material_grade` | string | - | Material grade, e.g. `"10.9"`. |
| `shank_length` | number | m | Length of the unthreaded shank. |
| `shank_diameter_ratio` | number or `"pitch"` | - | Shank diameter over nominal diameter. `"pitch"` sets the shank diameter equal to the thread basic pitch diameter (ISO 724), typical of reduced-shank studs. |
| `stud` | boolean | - | `true` for stud bolts, `false` for hexagonal-head bolts. |

**`MetricBolt`**: a bolt with custom properties (see `pyflange.bolts.MetricBolt`).

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `nominal_diameter` | number | m | Nominal thread diameter. |
| `thread_pitch` | number | m | Thread pitch. |
| `yield_stress` | number | Pa | Yield stress of the bolt material. |
| `ultimate_tensile_stress` | number | Pa | Ultimate tensile stress of the bolt material. |
| `elastic_modulus` | number | Pa | Young's modulus of the bolt material. |
| `poissons_ratio` | number | - | Poisson's ratio of the bolt material. |
| `shank_length` | number | m | Length of the unthreaded shank. |
| `shank_diameter_ratio` | number | - | Shank diameter over nominal diameter. |
| `stud` | boolean | - | `true` for stud bolts, `false` for hexagonal-head bolts. |

### Washers

**`ISOFlatWasher`**: ISO 7089 flat washer. Key: `designation` (e.g. `"M64"`).

**`FlatWasher`**: custom flat washer (see `pyflange.bolts.FlatWasher`).
Keys: `outer_diameter` (m), `inner_diameter` (m), `thickness` (m),
`elastic_modulus` (Pa), `poissons_ratio`.

### Nuts

**`ISOHexNut`**: ISO 4032 hexagonal nut. Key: `designation`.

**`RoundNut`**: standard flanged round nut. Key: `designation`.

**`HexNut`**: custom nut (see `pyflange.bolts.HexNut`). Keys:
`nominal_diameter` (m), `thickness` (m), `inscribed_diameter` (m),
`circumscribed_diameter` (m), `bearing_diameter` (m), `elastic_modulus` (Pa),
`poissons_ratio`.

### Gap

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `height` | number | m | Maximum gap height. |
| `angle` | number | rad | Angle subtended by the gap. |
| `shape_factor` | number | - | Gap shape factor (1.0 for a sinusoidal gap). Optional, default 1.0. |


BoltFatigueCurve
----------------

The S-N curve used in the deterministic analysis. The optional `type` key
selects the curve kind; if omitted, it is `"BoltFatigueCurve"`.

**`BoltFatigueCurve`**: IEC 61400-6 AMD1 bolt curve, with slopes 3 and 5,
knee point at 2·10⁶ cycles and size effect (see `pyflange.fatigue.BoltFatigueCurve`).

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `BoltNominalDiameter` | number | m | Bolt nominal diameter, for the size effect. |
| `DS_ref` | number | Pa | Reference stress range at 2·10⁶ cycles, before size effect and material factor. |
| `gamma_M` | number | - | Material safety factor. |

**`SingleSlopeFatigueCurve`**: keys `m`, `DS_ref` (Pa), `N_ref`.

**`DoubleSlopeFatigueCurve`**: keys `m1`, `m2`, `DS_ref` (Pa), `N_ref`,
where (`DS_ref`, `N_ref`) is the knee point.

In a Monte Carlo simulation, this section is not used: the random S-N curves
are defined by `Stats.FatigueClassFactor`.


MarkovMatrix
------------

The fatigue loads, as a list of bins of **flange (tower) bending moments**.

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `range` | list of numbers | N·m | Bending moment range of each bin. |
| `mean` | list of numbers | N·m | Mean bending moment of each bin. |
| `cycles` | list of numbers | - | Number of cycles of each bin. |
| `duration` | number | any | Duration of the load history (e.g. 30 for 30 years). The fatigue life is given in the same unit. |

The three lists must have the same length. Bins with zero range are allowed
and give no damage.


Stats
-----

The `Stats` section is required only for Monte Carlo simulations. It defines
the random variables used to generate random fatigue cases.

### Distributions

Random variables are described by *distribution objects*, selected by the
`Distribution` key:

| `Distribution` | Keys | Description |
|----------------|------|-------------|
| `"Normal"` | `Mean`, `CoV` | Normal distribution with given mean and coefficient of variation. |
| `"LogNormal"` | `Mean`, `CoV` | Log-normal distribution with given mean and coefficient of variation. |
| `"Weibull"` | `Shape`, `Scale` | Weibull distribution. |
| `"Beta"` | `Alpha`, `Beta`, `x0`, `x1` | Beta distribution defined over [`x0`, `x1`]. `x0` and `x1` are optional (default 0 and 1). |
| `"Product"` | `Dist1`, `Dist2` | Product of two independent variables, each described by one of the above distribution objects. Not allowed where a distribution is used as copula marginal. |

`Mean` is expressed in the unit of the variable.

### Stats keys

```json
"Stats": {
    "Preload":            { "Distribution": "Normal", "MeanRatio": 0.77, "CoV": 0.03 },
    "PlasticLoss":        { "Distribution": "Normal", "Mean": 0.03, "CoV": 0.33333 },
    "DFT":                { "Distribution": "Normal", "Mean": 250e-6, "CoV": 0.12 },
    "Gap":                { "JointDistribution": "IEC 61400-6:2020 AMD1", ... },
    "Tilt":               { "Distribution": "LogNormal", "Mean": 0.001745, "CoV": 0.5 },
    "LoadRangeFactor":    { "Distribution": "LogNormal", "Mean": 1.0, "CoV": 0.12 },
    "FatigueClassFactor": { "Distribution": "Normal", "Mean": 1.0, "CoV": 0.10, "DS_ref_mean": 62e6 },
    "AllowableDamage":    { "Distribution": "LogNormal", "Mean": 1.0, "CoV": 0.30 },
    "StressMultiplicationFactor": { "Distribution": "Product", "Dist1": { ... }, "Dist2": { ... } }
}
```

All keys are required.

| Key | Description |
|-----|-------------|
| `Preload` | Bolt preload before losses. `Distribution` is `"Normal"` or `"LogNormal"`, with `CoV` and a mean value given either as `MeanRatio` (fraction of `As·fy`, as `preload_ratio`) or as `Mean` (N, as `Fv`). |
| `PlasticLoss` | Preload loss due to plastic deformations, as a fraction of the preload (as `PPL`). Any distribution. |
| `DFT` | Coating dry film thickness, in m (as `DFT`). Any distribution. |
| `Gap` | Joint distribution of the gap parameters. See [Gap distributions](#gap-distributions). |
| `Tilt` | Flange tilt angle, in rad. Any distribution. Sampled values are reduced to the interval [0, π]. |
| `LoadRangeFactor` | Factor multiplying all the ranges of the Markov matrix. `"LogNormal"` only. |
| `FatigueClassFactor` | Factor multiplying the reference stress range of the bolt S-N curve. `"Normal"` only, with the additional key `DS_ref_mean` (Pa): the mean reference stress range. Each random curve is a `BoltFatigueCurve` with `DS_ref = DS_ref_mean × factor` and `gamma_M = 1.0`. |
| `AllowableDamage` | Allowable fatigue damage. Any distribution. |
| `StressMultiplicationFactor` | Stress multiplication factor applied to the bolt stress ranges. Any distribution, or a plain number for a constant factor. |

### Gap distributions

The `Gap` object is selected by the `JointDistribution` key.

**`"IEC 61400-6:2020 AMD1"`**: the gap angle is sampled from the `Angle`
distribution; the gap height is then sampled from the AMD1 log-normal
distribution for that angle (`pyflange.stats.gap_height_distribution`) and
capped at twice its 95% quantile (see `pyflange.stats.standard_gap_sampler`).

| Key | Description |
|-----|-------------|
| `FlatnessTolerance` | Flange flatness tolerance, in m/m (e.g. 0.0014 for 1.4 mm/m). |
| `Angle` | Distribution of the gap angle, in rad. Values outside (0, π) are discarded. |
| `ShapeFactor` | Distribution of the gap shape factor. |

**`"IEC 61400-6:2020 ED2"`**: the gap angle is sampled from the ED2 Beta
distribution and the gap height from the ED2 Weibull distribution
(see `pyflange.stats.standard_gap_sampler_ED2`). The shape factor is 1.0.

| Key | Description |
|-----|-------------|
| `FlatnessTolerance` | Flange flatness tolerance, in m/m. |

**`"ClaytonCopula"`**: gap height and gap angle are correlated via a Clayton
copula (see `pyflange.stats.clayton_gap_sampler`).

| Key | Description |
|-----|-------------|
| `Theta` | Clayton copula parameter (> 0). Kendall's rank correlation is θ/(θ+2). |
| `Height` | Marginal distribution of the gap height, in m (not `"Product"`). |
| `Angle` | Marginal distribution of the gap angle, in rad (not `"Product"`). |
| `ShapeFactor` | Distribution of the gap shape factor. |

### Monte Carlo model

Each random fatigue case is generated as follows (see
`pyflange.io.create_fatiguecase_sampler` and
`pyflange.stats.standard_PolynomialLFlangeSegment_sampler`):

- The flange segment has the same `type` and deterministic parameters as
  `FlangeSegment`, and random gap and tilt.
- The shell stiffness is always interpolated from FEA data for the random gap
  angle: the `FlangeSegment.k_shell` value is not used.
- The preload is sampled, reduced by the sampled plastic and settlement
  losses, and averaged over half of the bolts within the gap (at least one).
  The resulting segment has `PPL = 0` and `DFT = 0`, because the losses are
  already included in its preload.
- The Markov matrix ranges, the S-N curve, the allowable damage and the
  stress multiplication factor are sampled independently.

The `FlangeSegment.PPL`, `FlangeSegment.DFT`, `FlangeSegment.gap`,
`FlangeSegment.tilt_angle`, `FlangeSegment.k_shell`, `BoltFatigueCurve`,
`AllowableDamage` and `SMF` values are used in deterministic analyses only.


Result files
------------

The result files written by `pyflange determ` and `pyflange inspect` (and by
`pyflange.io.fatiguecase_to_dict`) use the same format, so they can be loaded
back as cases. They contain:

- `format_version`, `FlangeSegment`, `BoltFatigueCurve` and `MarkovMatrix`
  in the current format. The bolt, washer and nut are always written in
  their explicit forms (`MetricBolt`, `FlatWasher`, `HexNut`), the preload as
  `preload_ratio` and `k_shell` as a number. The `BoltFatigueCurve` section
  of a bolt curve also reports the design knee point (`m1`, `m2`, `DS12`,
  `N_ref`).
- `damage`, `allowable_damage`, `SMF`: the analysis results and settings.
- In verbose mode (always used by the CLI), `FatigueLife`, the bolt stress
  Markov matrix `BoltStressMarkovMatrix` (stresses in Pa) and a report of the
  flange segment intermediate quantities, added to the `FlangeSegment`
  section (`BoltPreload`, `Point1` to `Point4`, polynomial coefficients,
  `TransferFunctions`, stiffnesses, etc.).


Changes from version 0.1.0
--------------------------

Version 0.2.0 is a superset of version 0.1.0. The changes are:

- New `format_version` field (files without it are read as 0.1.0).
- `FlangeSegment.Fv` as an alternative to `preload_ratio`.
- `FlangeSegment.k_shell` accepts `"interp"`.
- Custom bolts (`MetricBolt`), washers (`FlatWasher`) and nuts (`HexNut`).
- `StandardMetricBolt.shank_diameter_ratio` accepts `"pitch"`.
- Optional `BoltFatigueCurve.type`, with `SingleSlopeFatigueCurve` and
  `DoubleSlopeFatigueCurve` curves.
- Optional top-level `AllowableDamage` and `SMF` for deterministic analyses.
- `Stats.Preload.Mean` as an alternative to `MeanRatio`.
- `Stats.StressMultiplicationFactor` accepts a plain number.
- **`Stats.FatigueClassFactor.DS_ref_mean` is required.** In version 0.1.0
  the mean reference stress range of the random S-N curves was implicitly
  62 MPa; version 0.1.0 files are upgraded with this value.
- `FlangeSegment.type` is optional, with default `"PolynomialLFlangeSegment"`.

> **Note on Monte Carlo simulations of ED2 cases.** The Forrester tool
> always generated random `PolynomialLFlangeSegment` (AMD1) segments, even
> for cases with `"type": "PolynomialLFlangeSegmentED2"`. PyFlange generates
> random segments of the declared type, so the Monte Carlo results of ED2
> cases differ from those computed with Forrester. Deterministic results
> and AMD1 cases are not affected.
