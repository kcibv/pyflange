# Validation

PyFlange is verified against the reference implementation of the polynomial
model by Siemens Gamesa (SGRE): a set of Excel workbooks that compute the
model of IEC 61400-6 AMD1, and of the coming ED2, step by step. The
comparison is automated and runs as part of the test suite, so every change
to the code is checked against the reference values.

## Validation cases

| Case | Model | Flange | Diameter | Gap | Tilt | Shell stiffness |
|------|-------|--------|----------|-----|------|-----------------|
| L1 | AMD1 | L, 120 × M80 | 7.5 m | sinusoidal | 0° | interpolated / simplified |
| L2 | AMD1 | L, 120 × M80 | 7.5 m | sinusoidal | 0.25° | interpolated / simplified |
| L3 | AMD1 | L, 120 × M80 | 7.5 m | shape factor 1.2 | 0° | interpolated / simplified |
| L4 | AMD1 | L | 4.6 m | sinusoidal | 0° | interpolated / simplified |
| T5 | AMD1 | T | 7.5 m | sinusoidal | 0° | interpolated / simplified |
| T6 | AMD1 | T | 7.5 m | shape factor 1.5 | 0° | interpolated / simplified |
| T7 | AMD1 | T | 7.5 m | shape factor 1.5 | 1° | interpolated / simplified |
| T8 | AMD1 | T | 7.5 m | sinusoidal | 1° | interpolated / simplified |
| L15 | ED2 | L, 192 × M64 | 8.0 m | sinusoidal | 0° | interpolated |

Each case is computed for the gap angles 30°, 60°, 90° and 120°. The AMD1
cases use the AMD1 design gap heights (95% quantile of the log-normal
distribution) and the classes `PolynomialLFlangeSegment` and
`PolynomialTFlangeSegment`. The ED2 case uses the ED2 design gap heights (75%
quantile of the Weibull distribution of `gap_height_distribution_ED2`), the
class `PolynomialLFlangeSegmentED2` and the S-N curve with material safety
factor γ<sub>M</sub> = 1.25. The 17 case files contain 3452 comparisons in
total: the bolt and flange stiffnesses, the gap height distribution, the
points and coefficients of the bolt force and moment polynomials, the gap
closing forces and the fatigue damage.

A quantity passes when it differs from the reference value by less than 1%.
For some quantities, each case file sets a larger tolerance, listed in its
`tolerances` object:

- the shell stiffness interpolated from the FEA database, and the gap
  stiffness and closing forces derived from it: up to 3%, mostly for the
  30° gap, because the workbooks interpolate the database differently;
- the fatigue damage: up to 5.5%, because it amplifies small differences of
  the bolt stress ranges;
- the bolt moment curve of T-flanges (initial slope `c_ini_M` and closed-gap
  moment `Ms_min`): up to 7%, because the reference workbooks fit a
  quadratic where PyFlange uses the exact cubic (see below);
- the initial slope of the bolt moment curve of the ED2 case (`c_ini_M`) for
  the 120° gap: 1.5%, because its reference value is small (about 25 times
  smaller than in the AMD1 cases);
- the slope of the bolt force polynomial at P1 (`c_ini_true`), a diagnostic
  quantity that is very sensitive to the polynomial coefficients: up to
  30%, and more in two cases where the reference value is close to zero.

Run the validation with `--failed` and a lower `--rtol` to see these
differences.

## Running the validation

The validation cases are in the `validations` folder of the
[source repository](https://github.com/kcibv/pyflange). After cloning it:

```doscon
py validations/validate.py validations/references/L1-interp.json
```

```
validations/references/L1-interp.json
    PASS  gap  30 deg  bolt.fy   [MPa]  expected      900   obtained      900   mismatch +0.000%   (allowable 1.00%)
    PASS  gap  30 deg  bolt.As   [mm²]  expected  4344.07   obtained  4344.05   mismatch -0.000%   (allowable 1.00%)
    ...
    PASS  gap  30 deg  Z2        [kN]   expected  1909.43   obtained  1920.17   mismatch +0.562%   (allowable 1.00%)
    ...
```

Without case files, all the cases are run: `py validations/validate.py`.

Useful options:

| Option | Effect |
|--------|--------|
| `--model ED2` | Run only the cases of a model version (`AMD1` or `ED2`). |
| `--flange-type T` | Run only the cases of a flange type (`L` or `T`). |
| `--list` | List the selected cases without running them. |
| `--failed` | Show only the comparisons that fail. |
| `--rtol 0.005` | Use a different allowable mismatch. |
| `--excel` | Also write the results into the PyFlange sheets of the SGRE workbooks of the cases, which then show the comparison in their `Comparison` sheet and the full list of checks in a `Validation` sheet (requires Excel and xlwings). The workbooks are updated in place; with `--excel DIR`, copies are saved in `DIR` instead. |

For example, `py validations/validate.py --model ED2 --excel` checks the ED2
case and writes its results to the ED2 workbook.

## Known differences

- The reference workbooks fit a quadratic polynomial through three points for
  the tensile bolt moment of T-flanges, while PyFlange evaluates the exact
  cubic of ref. [2], eq. G.61. The validation compares the quadratic
  coefficients of the workbooks with the quadratic through the PyFlange
  points.
- The ED2 workbooks still report the mean, coefficient of variation and 95%
  quantile of the AMD1 log-normal gap height distribution (k<sub>mean</sub>,
  COV<sub>k</sub>, k<sub>design</sub>), which are compared with the AMD1
  distribution of PyFlange, while the gap height used by the model is
  k<sub>used</sub>, the 75% quantile of the ED2 Weibull distribution.
- The ED2 workbooks report in the k<sub>fl</sub> row the flange stiffness
  divided by f<sub>tot</sub>, so that k<sub>gap,tot</sub> =
  f<sub>tot</sub> · (k<sub>shell,ini</sub> + k<sub>fl</sub>) gives the ED2 gap
  stiffness f<sub>tot</sub> · k<sub>shell,ini</sub> + k<sub>fl</sub>: the
  validation compares that row with k<sub>fl</sub> / f<sub>tot</sub>.
- The `PolynomialTFlangeSegmentED2` class is not covered by the reference
  workbooks yet. It is tested against the formulas of ref. [5] in the unit
  tests.

## Unit tests

Besides the validation, the test suite (`py -m pytest`) contains unit tests
of every module, and runs the example scripts of the
[tutorials](../tutorials/index.md).
