'''Export of the validation results to the SGRE validation workbooks.

For each validation case, the original workbook referenced by the case file is
copied to the output directory and the PyFlange results are written into its
PyFlange sheets (`PyFlange_Gap<angle>deg` in the L-flange workbooks,
`Gap<angle>deg` in the T-flange workbooks), exactly as the original scripts
`gen-cases.py` and `Case-D7500_T.py` did. Excel then recalculates the
comparison sheets (`Comparison` or `AMD1`) and the charts.

In addition, a `Validation` sheet is appended to each workbook, listing every
compared quantity with its mismatch, tolerance and pass/fail status, as
reported by `validate.py`.

The workbooks are edited through Excel (via xlwings), so that charts and form
controls are preserved: Excel must be installed.
'''

import os
import shutil
from math import inf

import numpy as np
import xlwings as xw

from pyflange.flangesegments import PolynomialLFlangeSegment


HERE = os.path.dirname(os.path.abspath(__file__))

mm = 0.001
kN = 1000
MPa = 1e6
GPa = 1e9
deg = np.pi/180


def pyflange_sheet_values (fseg, q):
    '''Values of the named cells of a PyFlange sheet.

    Args:
        fseg: the flange segment.
        q: the quantities evaluated by `validate.evaluate_quantities`.

    Returns:
        dict: maps the cell names to their values. Arrays are written in a row
            starting from the named cell. Names used only by the L-flange or
            only by the T-flange workbooks are both included: the names not
            present in a sheet are skipped.
    '''
    Z1 = fseg.shell_force_at_rest
    Z4 = fseg.shell_force_at_closed_gap
    Z = np.linspace(1.05*Z4, 1.05*fseg.shell_force_at_tensile_ULS, 100)

    values = {
        # Input values
        "a":                        fseg.a / mm,
        "b":                        fseg.b / mm,
        "shell_thickness":          fseg.s / mm,
        "t":                        fseg.t / mm,
        "central_angle":            fseg.central_angle / deg,
        "Radius":                   fseg.R / mm,
        "Z_dw":                     fseg.Zg / kN,
        "bolt.size":                fseg.bolt.designation,
        "bolt.pitch":               fseg.bolt.thread_pitch / mm,
        "bolt.Dn":                  fseg.bolt.nominal_diameter / mm,
        "bolt.Dsh":                 fseg.bolt.shank_diameter / mm,
        "bolt.Lsh":                 fseg.bolt.shank_length / mm,
        "bolt.fy":                  fseg.bolt.yield_stress / MPa,
        "bolt.fu":                  fseg.bolt.ultimate_tensile_stress / MPa,
        "bolt.E":                   fseg.bolt.elastic_modulus / GPa,
        "bolt.pretension":          fseg.Fv / kN,
        "Do":                       fseg.Do / mm,
        "Dw":                       (fseg.washer.outer_diameter if fseg.washer else fseg.nut.bearing_diameter) / mm,
        "gap.h":                    fseg.gap.height / mm,
        "gap.L":                    fseg.gap.angle * fseg.R / mm,
        "E_mod":                    fseg.E / GPa,
        "G_mod":                    fseg.G / GPa,

        # Polynomial data series
        "dataseries.Z":             list(Z / kN),
        "dataseries.Fs":            list(fseg.bolt_axial_force(Z) / kN),
        "dataseries.Ms":            list(fseg.bolt_bending_moment(Z)),

        # Polynomial data points
        "point1.Z":                 Z1 / kN,
        "point2.Z":                 fseg.shell_force_at_tensile_ULS / kN,
        "point3.Z":                 fseg.shell_force_at_small_displacement / kN,
        "point4.Z":                 Z4 / kN,
        "point1.Fs":                fseg.bolt_force_at_rest / kN,
        "point2.Fs":                fseg.bolt_force_at_tensile_ULS / kN,
        "point3.Fs":                fseg.bolt_force_at_small_displacement / kN,
        "point4.Fs":                fseg.bolt_axial_force(Z4) / kN,
        "point1.Ms":                fseg.bolt_moment_at_rest,
        "point2.Ms":                fseg.bolt_moment_at_tensile_ULS,
        "point3.Ms":                fseg.bolt_moment_at_small_displacement,
        "point4.Ms":                fseg.bolt_bending_moment(Z4),

        # Polynomial coefficients
        "Fs_coeff.tens.0":          list(fseg._tensile_force_polynomial.coef),
        "Fs_coeff.comp.0":          list(fseg._compressive_force_polynomial.coef),
        "Fs_coeff.const":           [fseg.bolt_axial_force(Z4), 0, 0],
        "Ms_coeff.tens.0":          list(fseg._tensile_moment_polynomial.coef),
        "Ms_coeff.comp.0":          list(fseg._compressive_moment_polynomial.coef),
        "Ms_coeff.const":           [fseg.bolt_bending_moment(Z4), 0, 0],

        # Miscellaneous model data
        "gap.angle":                fseg.gap.angle / deg,
        "bolt.As":                  fseg.bolt.thread_cross_section.area / mm**2,
        "a_prime":                  fseg.b / (fseg._prying_lever_ratio - 1) / mm,
        "Z0":                       fseg._ideal_shell_force_at_tensile_ULS / kN,
        "bolt.axial_stiffness":     fseg._bolt_axial_stiffness / (kN/mm),
        "bolt.bending_stiffness":   fseg._bolt_bending_stiffness / kN,
        "clamped_parts_stiffness":  fseg._flange_axial_stiffness / (kN/mm),
        "gap.stiffness":            fseg._gap_stiffness / 1e6,
        "DZ_gap_inclination":       fseg._tilt_neutralization_shell_force / kN,
        "DZ_gap":                   fseg._parallel_gap_neutralization_shell_force / kN,
        "DZ_gap_tot":               fseg._total_gap_neutralization_shell_force / kN,
        "stiffness_correction_factor": fseg._stiffness_correction_factor,
        "polynomial_initial_slope": fseg._polynomial_initial_slope,
        "true_force_initial_slope": q["c_ini_true"],
        "true_moment_initial_slope": q["c_ini_M"],

        # Model internal parameters
        "gap.k_mean":               q["gap.k_mean"],
        "gap.COV":                  q["gap.COV"],
        "gap.k_fac":                q["k_fac"],
        "gap.k_shell":              q["k_shell"],
        "gap.stiffness_correction_factor": q["f_tot"],     # L-flange workbooks
        "ftot":                     q["f_tot"],             # T-flange workbooks
        "flange.A_cf":              q["A_cf"],
        "flange.I_cf":              q["I_cf"],
        "gap.k_fl":                 q["k_fl"],
        "u":                        q["u"],
        "k_seg":                    q["k_seg"],
        "a_star":                   q["a_star"],
        "I_tg":                     q["I_tg"],
        "fatigue.damage":           q["damage"],            # L-flange workbooks
        "Damage":                   q["damage"],            # T-flange workbooks
    }

    if not isinstance(fseg, PolynomialLFlangeSegment):
        values["DZ_gap_c"] = fseg._early_prying_neutralization_shell_force / kN
        values["F_V_c"] = q["Fv_c"]

        # Bolt moment interpolation points at Z1 and Z1 + 100 N, as written by the
        # original script Case-D7500_T.py (see PolynomialTFlangeSegment._bolt_moment).
        s_avg = fseg.s * (1 + fseg.s_ratio) / 2
        c = fseg.central_angle * (fseg.R - s_avg/2)
        A_tg = c * fseg.t
        I_tg = c * fseg.t**3 / 12
        ak = fseg._stiffness_correction_factor
        Z0 = fseg._ideal_shell_force_at_tensile_ULS
        for i, Zi in enumerate([Z1, Z1 + 100], start=1):
            Fs = fseg.bolt_axial_force(Zi)
            phi_T_low = (Fs - fseg.Fv) * (fseg.b**2/(2*fseg.E*I_tg) + 1/(0.85*fseg.G*A_tg))
            phi_T_high = Zi * (fseg.b**2/(ak*4*fseg.E*I_tg) + 1/(0.85*fseg.G*A_tg))
            phi = phi_T_low + (phi_T_high - phi_T_low) / (Z0 - Z4) * (Zi - Z4)
            values[f"MZTlow{i}"] = phi_T_low * 2*fseg._bolt_bending_stiffness
            values[f"MZThigh{i}"] = phi_T_high * 2*fseg._bolt_bending_stiffness
            values[f"MZT{i}_"] = phi * 2*fseg._bolt_bending_stiffness
        values["DZdwDZ"] = (Z1 + 100) / kN

    return values


class ExcelExporter:
    '''Writes the validation results into copies of the SGRE validation workbooks.

    Args:
        out_dir (str): Directory where the workbooks are saved. Created if missing.
    '''

    def __init__ (self, out_dir):
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)
        self.app = xw.App(visible=False, add_book=False)
        self.app.display_alerts = False
        self.app.screen_updating = False

    def close (self):
        self.app.quit()

    def export (self, case, results):
        '''Exports the results of a case. Returns the path of the saved workbook.'''
        src = os.path.join(HERE, case["source"])
        dst = os.path.join(self.out_dir, os.path.basename(src))
        shutil.copyfile(src, dst)

        book = self.app.books.open(dst)
        try:
            for res in results:
                sheet = self._pyflange_sheet(book, res["angle"])
                self._write_named_cells(sheet, pyflange_sheet_values(res["fseg"], res["values"]))
            self._write_validation_sheet(book, case, results)
            self.app.calculate()
            book.save()
        finally:
            book.close()
        return dst

    @staticmethod
    def _pyflange_sheet (book, angle):
        names = [s.name for s in book.sheets]
        for name in (f"PyFlange_Gap{angle}deg", f"Gap{angle}deg"):
            if name in names:
                return book.sheets[name]
        raise KeyError(f"No PyFlange sheet for gap angle {angle} deg in {book.name}")

    @staticmethod
    def _write_named_cells (sheet, values):
        for name in sheet.names:
            key = name.name.split("!")[-1]
            if key not in values:
                continue
            try:
                cell = name.refers_to_range
            except Exception:
                continue    # broken name (#REF!)
            cell.value = values[key]

    @staticmethod
    def _write_validation_sheet (book, case, results):
        if "Validation" in [s.name for s in book.sheets]:
            book.sheets["Validation"].delete()
        sheet = book.sheets.add("Validation", after=book.sheets[-1])

        n_checks = sum(len(res["rows"]) for res in results)
        n_failed = sum(1 for res in results for r in res["rows"] if not r["passed"])
        sheet.range("A1").value = [
            [f"PyFlange validation - {case['id']}: {case['title']}"],
            [f"Result: {'PASS' if n_failed == 0 else 'FAIL'} - "
             f"{n_checks - n_failed}/{n_checks} checks within tolerance"],
            ["Mismatch = (PyFlange - SGRE) / |SGRE|"],
        ]
        sheet.range("A1").font.bold = True
        sheet.range("A1").font.size = 14

        header = ["Gap angle [deg]", "Quantity", "Symbol", "Unit", "SGRE", "PyFlange",
                  "Mismatch", "Tolerance", "Status"]
        table = [header]
        for res in results:
            for r in res["rows"]:
                mm_ = r["mismatch"] if r["mismatch"] not in (inf, -inf) else "inf"
                table.append([res["angle"], r["key"], r["symbol"], r["unit"], r["sgre"],
                              r["pyflange"], mm_, r["tolerance"], "PASS" if r["passed"] else "FAIL"])

        first_row = 5
        last_row = first_row + len(table) - 1
        sheet.range(f"A{first_row}").value = table
        sheet.range(f"A{first_row}:I{first_row}").font.bold = True
        sheet.range(f"A{first_row}:I{first_row}").color = (217, 225, 242)
        sheet.range(f"E{first_row+1}:F{last_row}").number_format = "0.000E+00"
        sheet.range(f"G{first_row+1}:H{last_row}").number_format = "0.000%"
        for i, row in enumerate(table[1:], start=first_row+1):
            if row[-1] == "FAIL":
                sheet.range(f"A{i}:I{i}").color = (255, 199, 206)
        sheet.range(f"A{first_row}:I{last_row}").api.AutoFilter(1)
        sheet.autofit("c")
