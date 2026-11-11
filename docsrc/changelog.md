# Changelog

The notable changes of each PyFlange version. The complete history is
available in the [GitHub repository](https://github.com/kcibv/pyflange/commits/main).

## 0.16.0 (unreleased)

!!! warning "Breaking changes"
    The ULS failure modes of the flange segments are now computed from
    properties of the segment instead of method arguments:

    - `failure_mode(fy_sh, fy_fl, gamma_0)` is replaced by the cached
      property `failure_mode` (the governing mode, a string) and by the
      properties `Zu_A`, `Zu_B`, `Zu_D`, `Zu_E` (the ultimate shell forces).
      Calling `fseg.failure_mode(...)` now raises
      `TypeError: 'str' object is not callable`.
    - `validate(fy_sh, fy_fl, gamma_0)` is now `validate()`, without
      arguments.
    - The yield stresses are no longer given directly: they are derived from
      the new field `steel_grade` (default `"S355"`) and from the thickness
      of shell and flange. The material factor is the new field `gamma_M0`
      (default 1.1).

    Migration:

    ```python
    # before
    mode, (Zu_A, Zu_B, Zu_D, Zu_E) = fseg.failure_mode(fy_sh=355e6, fy_fl=335e6)
    fseg.validate(fy_sh=355e6, fy_fl=335e6)

    # after (steel_grade="S355" and gamma_M0=1.1 are the defaults)
    fseg = PolynomialLFlangeSegment(..., steel_grade="S355", gamma_M0=1.1)
    mode = fseg.failure_mode
    Zu_A, Zu_B, Zu_D, Zu_E = fseg.Zu_A, fseg.Zu_B, fseg.Zu_D, fseg.Zu_E
    fseg.validate()
    ```

    To use other yield stresses than those of the table, set them before
    reading any ULS property: `fseg.fy_sh = 355e6; fseg.fy_fl = 335e6`.

- **ULS checks** (ref. [2], section G.1) for L- and T-flange segments,
  including the ED2 classes: cached properties `fy_sh`, `fy_fl`, `Zu_A`,
  `Zu_B`, `Zu_D`, `Zu_E`, `Zu` (the ultimate shell force F<sub>u</sub> of
  eq. G.1, the lowest of the four) and `failure_mode`.
- New fields `steel_grade` and `gamma_M0` of `PolynomialLFlangeSegment` and
  `PolynomialTFlangeSegment`, also available as optional keys of the JSON
  `FlangeSegment` section (format version unchanged) and as arguments of
  `stats.standard_PolynomialLFlangeSegment_sampler`.
- New `flangesegments.steel_yield_stress(grade, thickness)`: yield stress of
  S235, S275, S355, S460 and S500 for thicknesses up to 400 mm, after
  EN 10025-2.
- The verbose JSON report (`pyflange determ`, `pyflange inspect`) contains a
  new `ULS` object with the yield stresses, the ultimate shell forces and the
  governing failure mode.
- **ULS failure modes checked against IEC 61400-6:2020+AMD1:2025, Annex G.
  The following changes affect the ultimate shell forces:**
    - L-flanges, modes B, D and E: the plastic hinge M<sub>pl,3</sub> in the
      flange now uses the same lever arm as the hinge in the shell (b, or
      b'<sub>E</sub> for mode E), as in eqs. G.3 to G.5, instead of the
      shorter T-flange distance `b - s/2 - 0.8 r`. This lowers Zu by a few
      percent when the flange hinge governs (thick shell, thin flange); the
      previous values were non-conservative.
    - Modes B, D and E are limited to the plastic resistance of the shell in
      tension and of the flange in shear, above which eq. G.11 is not
      defined. Before, the ultimate force could exceed them (e.g. the L-flange
      of the test suite had Zu_E above the shell tensile resistance).
    - The washer diameter D of eqs. G.9 and G.29 is limited to "width
      across flats of the nut + 2 × washer height" for narrow washers
      (Figure G.3).
    - The shell resistance uses the thinner of the top and bottom shell
      (`s_ratio < 1`), as G.2.1 a).
    - The fixed-point equations are solved by bisection instead of `fsolve`,
      which is robust to the kink of the interaction formula.
- The code shared by `PolynomialLFlangeSegment` and
  `PolynomialTFlangeSegment` has been moved to `PolynomialFlangeSegment`,
  with no change in the results (all the quantities of the L-, T- and ED2
  segments were verified to be bit-identical before and after). This only
  matters for code that subclasses `PolynomialFlangeSegment` directly: it is
  now a dataclass with the fields common to L- and T-flanges, it implements
  the reference points of the polynomials (ref. [2], Table G.1), and the
  child classes must implement the abstract quantities listed in its
  docstring.
- The API pages show the members inherited from the base classes (e.g.
  `bolt_axial_force`, `Zu` and `validate` on the L- and T-flange pages).

## 0.15.0 (2026-11-02)

- `PolynomialLFlangeSegmentED2` validated against the SGRE ED2 reference
  workbook (case L15, D8000 L-flange with ED2 design gap heights).
- Validation script: without case files it runs all the cases; new options
  `--model`, `--flange-type` and `--list`. The L-flange segments are built
  through `pyflange.io`, so the validation covers the JSON case code path.
- **Fixed: interpolated shell stiffness for `s_ratio != 1`.** With
  `"k_shell": "interp"` in a JSON case (`create_flangesegment`, the Monte
  Carlo sampler, `create_realization`, `pyflange determ --gap-angle`) and in
  `standard_PolynomialLFlangeSegment_sampler`, the shell stiffness is now
  interpolated for the average thickness of the top and bottom shell,
  `s * (1 + s_ratio) / 2`, as in the simplified formula (eqs. G.43/G.44)
  and in the SGRE validation workbooks, instead of `s`. **This changes the
  results whenever `s_ratio != 1`**: k<sub>shell</sub> is higher and the
  damages differ (e.g. D 8 m, s = 37.2 mm, s_ratio = 2.15: k<sub>shell</sub>
  at 30° from 1853 to 2942 MN/m, ED2 damage at 90°–120° up to 28 % higher).
  Results with `s_ratio = 1` are unchanged.
- New `flangesegments.average_shell_thickness`, the thickness to pass to
  `shell_stiffness`.
- `stats.standard_markov_matrix_sampler`: the default CoV of the load range
  factor X<sub>FLS</sub> is now 0.168 instead of 0.12 (refs. [1], Table 7,
  and [5], Table 8). Calls passing `range_CoV` explicitly and JSON cases
  (which always give the `LoadRangeFactor` CoV) are unchanged. The example
  case `fatigue-case.json` and the Monte Carlo tutorials now use 0.168 too.

## 0.14.0 (2026-10-23)

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
- Removed `utils.log_data` and `utils.read_data_log`: the intermediate
  quantities of the flange segment models (k<sub>fac</sub>, f<sub>tot</sub>,
  A<sub>cf</sub>, I<sub>cf</sub>, u, k<sub>seg</sub>, I<sub>tg</sub>,
  F<sub>V,c</sub>, ...) are now private cached properties, computed on
  demand. The verbose JSON output no longer contains the `gap_k_mean` and
  `gap_COV` keys, which were always `null`.

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
