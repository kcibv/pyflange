'''Validation of PyFlange against the SGRE reference values.

Each validation case is defined by a JSON file in `validations/references/`
(generated from the SGRE validation workbooks by `extract_references.py`). For
every gap angle of the case, `validate` builds the corresponding PyFlange
flange segment, evaluates the quantities listed in the case file, compares
them with the SGRE reference values and yields the result of each comparison.

A quantity passes when its relative mismatch `(pyflange - sgre) / |sgre|` does
not exceed the allowable mismatch in absolute value. The allowable mismatch is
`rtol` (default 1%), unless the case file overrides it in its optional
`tolerances` object, which maps quantity keys (or `"default"`) either to a
relative tolerance or to a list of relative tolerances, one per gap angle.

From the command line, the script displays the comparisons of the given case
files and exits with code 0 if all of them pass, 1 otherwise. With `--excel`,
the results are also exported to copies of the original SGRE validation
workbooks (see `excel_export.py`, requires xlwings and Excel).

Usage:
    py validations/validate.py PATH [PATH ...] [--rtol RTOL] [--failed] [--excel [DIR]]

Examples:
    py validations/validate.py validations/references/L1-interp.json
    py validations/validate.py references/L1-interp.json     # also relative to validations/
    py validations/validate.py validations/references/*.json --failed
    py validations/validate.py validations/references/T5-simpl.json --rtol 0.005
    py validations/validate.py validations/references/*.json --excel
'''

import os
import sys
import json
import logging
import argparse
from math import pi, inf

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EXCEL_DIR = os.path.join(HERE, "results")

# Make sure the pyflange in this repository is the one being validated
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "src"))

from pyflange.flangesegments import (PolynomialLFlangeSegment, PolynomialTFlangeSegment,
                                     Gap, shell_stiffness, bolt_markov_matrix)
from pyflange.bolts import MetricBolt, HexNut, FlatWasher
from pyflange.fatigue import MarkovMatrix, BoltFatigueCurve, BoltFatigueAnalysis
from pyflange.stats import gap_height_distribution


# Units of measurement
mm = 0.001
kN = 1000
MPa = 1e6
deg = pi/180

DEFAULT_RTOL = 0.01



# ------------------------------------------------------------------------------
#  Model construction
# ------------------------------------------------------------------------------

_markov_cache = {}

def load_markov_matrix (path):
    path = os.path.abspath(path)
    if path not in _markov_cache:
        with open(path) as f:
            data = json.load(f)
        _markov_cache[path] = MarkovMatrix(
            range = np.array(data["range"]),
            mean = np.array(data["mean"]),
            cycles = np.array(data["cycles"]),
            duration = data["duration"])
    return _markov_cache[path]


def create_flange_segment (case, gap_angle):
    '''Creates the flange segment of a validation case, for the given gap angle (rad).

    Returns the flange segment and the gap height distribution.
    '''
    seg = case["segment"]
    R = seg["D"] / 2

    gap_dist = gap_height_distribution(seg["D"], case["gap"]["flatness_tolerance"], gap_angle * R)
    gap = Gap(height = gap_dist.ppf(case["gap"]["fractile"]),
              angle = gap_angle,
              shape_factor = case["gap"]["shape_factor"])

    if case["shell_stiffness"] == "interpolated":
        k_shell = float(shell_stiffness(R, seg["s"] * (1 + seg["s_ratio"]) / 2, gap_angle))
    elif case["shell_stiffness"] == "simplified":
        k_shell = None
    else:
        raise ValueError(f"Unknown shell stiffness option '{case['shell_stiffness']}'")

    FlangeSegment = {"L": PolynomialLFlangeSegment, "T": PolynomialTFlangeSegment}[seg["flange_type"]]
    fseg = FlangeSegment(
        a = seg["a"],
        b = seg["b"],
        s = seg["s"],
        t = seg["t"],
        R = R,
        central_angle = 2*pi / seg["n_bolts"],
        Zg = seg["Zg"],
        bolt = MetricBolt(**seg["bolt"]),
        Fv = seg["Fv"],
        Do = seg["Do"],
        washer = FlatWasher(**seg["washer"]) if seg["washer"] else None,
        nut = HexNut(**seg["nut"]),
        gap = gap,
        tilt_angle = case["tilt_angle_deg"] * deg,
        s_ratio = seg["s_ratio"],
        k_shell = k_shell)

    return fseg, gap_dist


def fatigue_damage (fseg, case):
    '''Bolt fatigue damage according to the fatigue parameters of the case.'''
    fat = case["fatigue"]
    # The Markov matrix path is relative to the case file
    flange_mkvm = load_markov_matrix(os.path.join(case["_dir"], fat["markov_matrix"]))

    # Markov matrix bins with zero range give a harmless division by zero in the S-N curve
    with np.errstate(divide="ignore"):
        if fat["bending_factor"] == "IEC":
            # Standard IEC 61400-6 AMD1 analysis
            if fat["mean_factor"] != 1.0 or fat["range_factor"] != 1.0:
                raise ValueError("Mean and range factors are not supported with the 'IEC' bending factor")
            return BoltFatigueAnalysis(fseg, flange_mkvm).damage

        bolt_mkvm = bolt_markov_matrix(fseg, flange_mkvm,
                                       bending_factor = fat["bending_factor"],
                                       mean_factor = fat["mean_factor"],
                                       range_factor = fat["range_factor"])
        return BoltFatigueCurve(fseg.bolt.nominal_diameter).cumulated_damage(bolt_mkvm)



# ------------------------------------------------------------------------------
#  Compared quantities
#
#  All values are expressed in the units of the SGRE workbooks.
# ------------------------------------------------------------------------------

def evaluate_quantities (fseg, gap_dist, case):
    '''Returns a dict with the value of each quantity known to the validation.'''
    is_L = isinstance(fseg, PolynomialLFlangeSegment)

    Z1 = fseg.shell_force_at_rest
    Z4 = fseg.shell_force_at_closed_gap
    Ks = fseg._bolt_axial_stiffness
    Kp = fseg._flange_axial_stiffness
    Fs_poly = fseg._tensile_force_polynomial
    Ms_poly = fseg._tensile_moment_polynomial

    M1 = fseg.bolt_moment_at_rest
    M2 = fseg.bolt_moment_at_tensile_ULS
    M3 = fseg.bolt_moment_at_small_displacement
    Ms_min = fseg.bolt_bending_moment(Z4)

    q = {
        "bolt.fy":              fseg.bolt.yield_stress / MPa,
        "bolt.As":              fseg.bolt.thread_cross_section.area / mm**2,
        "bolt.Fv":              fseg.Fv / kN,
        "bolt.dS":              1 / Ks,
        "flange.dP":            1 / Kp,
        "p":                    Ks / (Ks + Kp),
        "gap.angle":            fseg.gap.angle / deg,
        "gap.L":                fseg.gap.angle * fseg.R,
        "D":                    2 * fseg.R,
        "gap.k_mean":           gap_dist.mean() / mm,
        "gap.COV":              gap_dist.std() / gap_dist.mean(),
        "gap.height":           fseg.gap.height / mm,
        "Z1":                   Z1 / kN,
        "Fs1":                  fseg.bolt_force_at_rest / kN,
        "Z2":                   fseg.shell_force_at_tensile_ULS / kN,
        "Fs2":                  fseg.bolt_force_at_tensile_ULS / kN,
        "Z3":                   fseg.shell_force_at_small_displacement / kN,
        "Fs3":                  fseg.bolt_force_at_small_displacement / kN,
        "Fs_poly.a2":           Fs_poly.coef[2] * kN,
        "Fs_poly.a1":           Fs_poly.coef[1],
        "Fs_poly.a0":           Fs_poly.coef[0] / kN,
        "DZ_gap":               -fseg._parallel_gap_neutralization_shell_force / kN,
        "I_tg":                 fseg._I_tg / mm**4,
        "A_cf":                 fseg._A_cf / mm**2,
        "I_cf":                 fseg._I_cf / mm**4,
        "k_gap_tot":            fseg._gap_stiffness / 1e6,
        "f_tot":                fseg._f_tot,
        "k_shell":              fseg._k_shell / 1e6,
        "k_fac":                fseg._k_fac,
        "k_fl":                 fseg._k_flange / 1e6,
        "DZ_gap_inclination":   fseg._tilt_neutralization_shell_force / kN,
        "DZ_gap_tot":           -fseg._total_gap_neutralization_shell_force / kN,
        "Z0":                   fseg._ideal_shell_force_at_tensile_ULS / kN,
        "ak":                   fseg._stiffness_correction_factor,
        "k_seg":                fseg._k_seg / 1e6,
        "u":                    fseg._u / mm,
        "Z_dw":                 fseg.Zg / kN,
        "Z_close":              Z4 / kN,
        "c_ini_true":           Fs_poly.deriv()(Z1),
        "Fs_min":               fseg.bolt_axial_force(Z4) / kN,
        "c_ini_mod":            fseg._polynomial_initial_slope,
        "a_star":               fseg._a_star / mm,
        "bolt.bS":              1 / (fseg._bolt_bending_stiffness * 1000),
        "c_ini_M":              Ms_poly.deriv()(Z1) * kN,
        "Ms_min":               Ms_min,
        "M1":                   M1,
        "M2":                   M2,
        "M3":                   M3,
        "Ms_poly.a2":           Ms_poly.coef[2] * kN**2,
        "Ms_poly.a1":           Ms_poly.coef[1] * kN,
        "Ms_poly.a0":           Ms_poly.coef[0],
        "damage":               fatigue_damage(fseg, case),
    }

    if not is_L:
        # PyFlange evaluates the T-flange bolt moment exactly, as the cubic polynomial
        # of IEC 61400-6, eq. G.61, while the SGRE sheets fit a quadratic polynomial
        # through the points P1, P2 and P3. The SGRE coefficients are therefore compared
        # with the quadratic polynomial through the same points of PyFlange.
        Z = [Z1, fseg.shell_force_at_tensile_ULS, fseg.shell_force_at_small_displacement]
        Ms_fit = np.polynomial.polynomial.polyfit(Z, [M1, M2, M3], 2)
        q["Ms_poly.a2"] = Ms_fit[2] * kN**2
        q["Ms_poly.a1"] = Ms_fit[1] * kN
        q["Ms_poly.a0"] = Ms_fit[0]

        # The SGRE T-flange sheets report the tilt contribution as a positive
        # reduction, like DZgap and DZgap,tot (the L-flange sheets keep the sign).
        q["DZ_gap_inclination"] = -q["DZ_gap_inclination"]
        q["DZ_gap_c"] = -fseg._early_prying_neutralization_shell_force / kN
        q["Fv_c"] = fseg._Fv_c / kN
        q["Ms_T.at_rest"] = fseg._bolt_moment(Z1, fseg.bolt_axial_force(Z1))
        q["Ms_T.at_rest+100N"] = fseg._bolt_moment(Z1 + 100, fseg.bolt_axial_force(Z1 + 100))

    return q



# ------------------------------------------------------------------------------
#  Comparison and reporting
# ------------------------------------------------------------------------------

def mismatch (value, reference):
    '''Relative mismatch between a PyFlange value and its reference.'''
    if reference == 0:
        return 0.0 if abs(value) < 1e-9 else inf
    return (value - reference) / abs(reference)


def tolerance (case, key, angle_index, rtol):
    '''Relative tolerance of a quantity at a gap angle (case overrides, else `rtol`).'''
    tolerances = case.get("tolerances", {})
    tol = tolerances.get(key, tolerances.get("default", rtol))
    return tol[angle_index] if isinstance(tol, list) else tol


def load_case (path):
    '''Loads a validation case file.

    Relative paths are looked up in the current directory first, then in the
    `validations/` directory, so that e.g. `references/L1-interp.json` works
    from both.
    '''
    if not os.path.isabs(path) and not os.path.exists(path):
        path = os.path.join(HERE, path)
    with open(path, encoding="utf-8") as f:
        case = json.load(f)
    case["_dir"] = os.path.dirname(os.path.abspath(path))
    return case


def run_case (case, rtol=DEFAULT_RTOL):
    '''Runs a validation case.

    Returns a list with one item per gap angle: a dict with the gap angle, the
    flange segment, the gap height distribution, the evaluated quantities and
    the comparison rows (one dict per compared quantity).
    '''
    results = []
    for i, angle in enumerate(case["gap_angles_deg"]):
        fseg, gap_dist = create_flange_segment(case, angle * deg)
        values = evaluate_quantities(fseg, gap_dist, case)

        rows = []
        for key, ref in case["references"].items():
            tol = tolerance(case, key, i, rtol)
            mm_ = mismatch(values[key], ref["values"][i])
            rows.append(dict(key=key, symbol=ref["symbol"], unit=ref["unit"],
                             sgre=ref["values"][i], pyflange=values[key],
                             mismatch=mm_, tolerance=tol, passed=abs(mm_) <= tol))

        results.append(dict(angle=angle, fseg=fseg, gap_dist=gap_dist, values=values, rows=rows))
    return results


def validate (path, rtol=DEFAULT_RTOL):
    '''Validates PyFlange against the reference values of a validation case.

    Args:
        path (str): Path of the validation case JSON file (e.g. `references/L1-interp.json`).
        rtol (float): Allowable relative mismatch of the quantities that have no
            tolerance in the case file. Default: 1%.

    Yields:
        dict: One comparison per reference quantity and gap angle, with keys:

            - `case`: ID of the validation case;
            - `gap_angle`: gap angle in degrees;
            - `parameter`: key of the compared quantity (e.g. `"Fs_poly.a2"`);
            - `symbol`, `unit`: symbol and unit of the quantity in the SGRE workbooks;
            - `expected`: SGRE reference value;
            - `obtained`: PyFlange value;
            - `mismatch`: relative mismatch `(obtained - expected) / |expected|`;
            - `allowable_mismatch`: maximum allowed absolute value of the mismatch;
            - `pass`: `True` if `|mismatch| <= allowable_mismatch`, `False` otherwise.
    '''
    case = load_case(path)
    for res in run_case(case, rtol):
        for r in res["rows"]:
            yield {
                "case":                 case["id"],
                "gap_angle":            res["angle"],
                "parameter":            r["key"],
                "symbol":               r["symbol"],
                "unit":                 r["unit"],
                "expected":             r["sgre"],
                "obtained":             r["pyflange"],
                "mismatch":             r["mismatch"],
                "allowable_mismatch":   r["tolerance"],
                "pass":                 r["passed"],
            }


def fmt_pct (x):
    return "inf" if x in (inf, -inf) else f"{x:+.3%}"


def print_comparison (c):
    flag = "PASS" if c["pass"] else "FAIL"
    print(f"    {flag}  gap {c['gap_angle']:>3d} deg  {c['parameter']:<20s} {'[' + c['unit'] + ']':<11s}"
          f" expected {c['expected']:>14.6g}   obtained {c['obtained']:>14.6g}"
          f"   mismatch {fmt_pct(c['mismatch']):>9s}   (allowable {c['allowable_mismatch']:.2%})")



# ------------------------------------------------------------------------------
#  Command line interface
# ------------------------------------------------------------------------------

def main (argv=None):
    parser = argparse.ArgumentParser(
        description="Validate PyFlange against the SGRE reference values of validation cases.")
    parser.add_argument("paths", nargs="+", metavar="PATH",
                        help="validation case JSON file(s), e.g. 'references/L1-interp.json'")
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL,
                        help=f"allowable relative mismatch of the quantities without a tolerance "
                             f"in the case file (default: {DEFAULT_RTOL})")
    parser.add_argument("--failed", action="store_true",
                        help="display only the failed comparisons")
    parser.add_argument("--excel", nargs="?", const=DEFAULT_EXCEL_DIR, metavar="DIR",
                        help="also export the results to copies of the SGRE validation workbooks "
                             f"in DIR (default: {os.path.relpath(DEFAULT_EXCEL_DIR)}). Requires Excel.")
    parser.add_argument("--show-warnings", action="store_true",
                        help="show the warnings logged by pyflange")
    args = parser.parse_args(argv)

    if not args.show_warnings:
        logging.disable(logging.WARNING)

    summary = []
    for path in args.paths:
        print(f"\n{path}")
        n_total = n_failed = 0
        for c in validate(path, args.rtol):
            n_total += 1
            n_failed += not c["pass"]
            if not args.failed or not c["pass"]:
                print_comparison(c)
        print(f"    {n_total - n_failed}/{n_total} comparisons passed")
        summary.append((path, n_failed))

    if args.excel:
        from excel_export import ExcelExporter
        exporter = ExcelExporter(args.excel)
        try:
            for path in args.paths:
                case = load_case(path)
                print(f"{case['id']} exported to {exporter.export(case, run_case(case, args.rtol))}")
        finally:
            exporter.close()

    failed = [(path, n) for path, n in summary if n > 0]
    print("\n" + "=" * 80)
    print(f"{len(summary) - len(failed)}/{len(summary)} cases passed.")
    for path, n in failed:
        print(f"    {path}: {n} failed comparisons")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
