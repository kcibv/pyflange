# Validation

PyFlange is verified against the reference implementation of the polynomial
model by Siemens Gamesa (SGRE): a set of Excel workbooks that compute the
model of IEC 61400-6 AMD1 step by step. The comparison is automated and runs
as part of the test suite, so every change to the code is checked against
the reference values.

## Validation cases

| Case | Flange | Diameter | Gap | Tilt | Shell stiffness |
|------|--------|----------|-----|------|-----------------|
| L1 | L, 120 × M80 | 7.5 m | sinusoidal | 0° | interpolated / simplified |
| L2 | L, 120 × M80 | 7.5 m | sinusoidal | 0.25° | interpolated / simplified |
| L3 | L, 120 × M80 | 7.5 m | shape factor 1.2 | 0° | interpolated / simplified |
| L4 | L | 4.6 m | sinusoidal | 0° | interpolated / simplified |
| T5 | T | 7.5 m | sinusoidal | 0° | interpolated / simplified |
| T6 | T | 7.5 m | shape factor 1.5 | 0° | interpolated / simplified |
| T7 | T | 7.5 m | shape factor 1.5 | 1° | interpolated / simplified |
| T8 | T | 7.5 m | sinusoidal | 1° | interpolated / simplified |

Each case is computed for the gap angles 30°, 60°, 90° and 120°, with the
AMD1 design gap heights. The 16 case files contain 3264 comparisons in total:
the bolt and flange stiffnesses, the gap height distribution, the points and
coefficients of the bolt force and moment polynomials, the gap closing forces
and the fatigue damage.

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

Useful options:

| Option | Effect |
|--------|--------|
| `--failed` | Show only the comparisons that fail. |
| `--rtol 0.005` | Use a different allowable mismatch. |
| `--excel` | Also write the results to copies of the SGRE workbooks in `validations/results/` (requires Excel and xlwings). |

To run all the cases: `py validations/validate.py validations/references/*.json`.

## Known differences

- The reference workbooks fit a quadratic polynomial through three points for
  the tensile bolt moment of T-flanges, while PyFlange evaluates the exact
  cubic of ref. [2], eq. G.61. The validation compares the quadratic
  coefficients of the workbooks with the quadratic through the PyFlange
  points.
- The ED2 classes are not covered by the reference workbooks. They are
  tested against the formulas of ref. [5] in the unit tests.

## Unit tests

Besides the validation, the test suite (`py -m pytest`) contains unit tests
of every module, and runs the example scripts of the
[tutorials](../tutorials/index.md).
