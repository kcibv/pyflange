# Changelog

The notable changes of each PyFlange version. The complete history is
available in the [GitHub repository](https://github.com/kcibv/pyflange/commits/main).

## 0.13.0 (unreleased)

- **Model aligned with the published IEC 61400-6:2020+AMD1:2025**: docstrings
  cite the published equation numbers; T-flange ULS failure modes D and E,
  initial slope bound and b'<sub>E</sub> corrected; exact cubic bolt moment
  for T-flanges (eq. G.61); Tobinaga & Ishihara correction applied only
  within its validity range.
- T-flanges: preload losses (`PPL`, `DFT`, `settlement_factor`), callable
  `k_shell`, and `tilt_angle` is now a required argument.
- New `PolynomialTFlangeSegmentED2`.
- New **JSON case format** (version 0.2.0) and `pyflange.io` module.
- New **command-line interface**: `pyflange determ`, `sample`, `inspect`,
  `stats`.
- ED2 gap statistics: `gap_height_distribution_ED2`,
  `gap_angle_distribution_ED2`, `standard_gap_sampler_ED2`.
- New samplers: `product_sampler`, `clayton_copula_sampler`,
  `clayton_gap_sampler`.
- `standard_hole_diameter`; `StandardMetricBolt` accepts
  `shank_diameter_ratio="pitch"`.
- Automated validation against the SGRE reference workbooks.
- New documentation website.

## 0.12.0 (2026-10-09)

- Documentation of `README.md`, `bolts` and `fatigue` improved.

## 0.11.4 (2026-05-13)

- `standard_gap_sampler`: gap angles larger than 180° are discarded and
  drawn again.

## 0.11.3 (2026-05-04)

- `standard_gap_sampler`: realizations with a gap height above the cap are
  discarded and drawn again, instead of being set to the 95% quantile.

## 0.11.2 (2026-04-16)

- First draft of the ED2 L-flange model (`PolynomialLFlangeSegmentED2`).
- New `weibull_sampler` and `beta_sampler`.

## 0.11.1 (2026-03-16)

- Gap height cap of `standard_gap_sampler` as a parameter.
- Stress multiplication factor in `fatigue_case_sampler`, constant or random.
- Fixed `reduce_to_pi`.

## 0.11.0 (2026-03-13)

- Stress multiplication factor (`SMF`) for the bolt stress ranges.
- `standard_bolt_fatigue_curve_sampler`: mean reference stress range as a
  parameter.
- Fixed the averaging of the bolt preloads in the flange segment sampler.
- Removed `markov_matrix_from_SGRE_format`.

## 0.10.0 (2026-03-06)

- Bolt preload losses due to settlement and plastic strain (`PPL`, `DFT`,
  `settlement_factor`) in `PolynomialLFlangeSegment` and in the flange
  segment sampler.

## 0.9.1 (2025-12-19)

- Minor fixes and API documentation updates.

## 0.9.0 (2025-12-16)

- New `Gap` class: the flange segment parameters `gap_height`, `gap_angle`
  and `gap_shape_factor` are replaced by a single `gap` parameter.

## 0.8.0 (2025-11-12)

- New `stats` module, with probability distributions and random samplers for
  Monte Carlo simulations; `gap` module deprecated.
- New `BoltFatigueAnalysis` class.
- `k_shell` can be a function of the flange segment.

## 0.7.4 (2025-10-21)

- New `MarkovMatrix` class, with an explicit load `duration`.
- Failure mode D uses c' instead of c.

## 0.7.3 (2025-09-10)

- Fixed the iterative calculation of failure modes D and E.
- Fixed the initial slope of the L-flange bolt force curve.

## 0.7.2 (2025-08-20)

- `PolynomialTFlangeSegment.validate` requires failure mode A instead of B.

## 0.7.0 - 0.7.1 (2025-05)

- Documentation prepared for the public release.
- Bolt, nut and washer databases moved to CSV files.

## 0.6.0 (2024-11-13)

- Gap stiffness factor as a function of the gap angle.
- Shell stiffness interpolated from FEA results (`shell_stiffness`).
- `bolt_markov_matrix` function; `BoltCrossSection` class.

## 0.5.0 (2024-09-16)

- New `fatigue` module with `BoltFatigueCurve`.
- Washers and nuts in the flange segment models.

## 0.4.0 (2024-06-26)

- `PolynomialTFlangeSegment`: gap shapes other than sinusoidal.

## 0.3.0 (2024-06-18)

- `PolynomialLFlangeSegment`: flange tilt, validated against the SGRE
  reference implementation.

## 0.1.0 - 0.2.0 (2024)

- First versions: polynomial L-flange segment model and bolts.
