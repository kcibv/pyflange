
from pyflange.flangesegments import PolynomialLFlangeSegment
from pyflange.bolts import MetricBolt
from pyflange.gap import gap_height_distribution

from math import pi
import numpy as np

import os
from xlwings import Book

# Units of measurement
m = 1
mm = 0.001*m

kg = 1
t = 1000*kg

s = 1

N = kg*m/s**2
kN = 1000*N

Pa = 1
MPa = 1e6*Pa
GPa = 1e9*Pa

rad = 1
deg = (pi/180)*rad



# Bolt

M80 = MetricBolt(
    nominal_diameter = 80*mm,
    thread_pitch = 6*mm,
    shank_diameter_ratio = 76.1/80,
    shank_length = 270*mm,
    yield_stress = 900*MPa,
    ultimate_tensile_stress = 1000*MPa,
    stud = True)




# Polinomial Segment Model

def create_flange_segment (gap_angle):

    D = 7500*mm
    t_sh = 72*mm
    n = 120 # number of bolts
    gap_length = gap_angle * D/2
    gap = gap_height_distribution(D, 0.0014, gap_length)

    k_mean = gap.mean()
    COV_k = gap.std() / k_mean

    print(f"gap mean height: k_mean = {k_mean/mm:.3f} mm")
    print(f"gap COV: COV_k = {COV_k:.3f}")

    fseg = PolynomialLFlangeSegment(

        a = 232.5*mm,           # distance between inner face of the flange and center of the bolt hole
        b = 166.5*mm,           # distance between center of the bolt hole and center-line of the shell
        s = t_sh,           # shell thickness
        t = 200.0*mm,           # flange thickness
        c = 2*pi/n * (D - t_sh)/2,   # shell arc length
        R = D/2,          # shell outer curvature radius

        Zg = -14795*kN / n, # load applied to the flange segment shell at rest
                                                # (normally dead weight of tower + RNA, divided by the number of bolts)

        bolt = M80,
        Fv = 2876*kN,                            # applied bolt preload

        Do = 86*mm,     # bolt hole diameter
        Dw = 140*mm,    # washer diameter

        gap_height = gap.ppf(0.95),   # maximum longitudinal gap height
        gap_angle = gap_angle,  # longitudinal gap length

        s_ratio = 102/72        # ratio of bottom shell thickness over tower shell thickness
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    return fseg




# Write Results to Excel
def set_cell_value (book, name, value):
    book.names[name].refers_to_range.value = value


def flange_segment_model_to_excel (book, sheet_name, fseg):

    # Write input values to excel
    set_cell_value(book, f"{sheet_name}!a", fseg.a/mm)
    set_cell_value(book, f"{sheet_name}!b", fseg.b/mm)
    set_cell_value(book, f"{sheet_name}!shell_thickness", fseg.s/mm)
    set_cell_value(book, f"{sheet_name}!t", fseg.t/mm)
    set_cell_value(book, f"{sheet_name}!shell_arc_length", fseg.c/mm)
    set_cell_value(book, f"{sheet_name}!Radius", fseg.R/mm)
    set_cell_value(book, f"{sheet_name}!Z_dw", fseg.Zg/kN)
    set_cell_value(book, f"{sheet_name}!bolt.size", fseg.bolt.designation)
    set_cell_value(book, f"{sheet_name}!bolt.pitch", fseg.bolt.thread_pitch/mm)
    set_cell_value(book, f"{sheet_name}!bolt.Dn", fseg.bolt.nominal_diameter/mm)
    set_cell_value(book, f"{sheet_name}!bolt.Dsh", fseg.bolt.shank_diameter/mm)
    set_cell_value(book, f"{sheet_name}!bolt.Lsh", fseg.bolt.shank_length)
    set_cell_value(book, f"{sheet_name}!bolt.fy", fseg.bolt.yield_stress/MPa)
    set_cell_value(book, f"{sheet_name}!bolt.fu", fseg.bolt.ultimate_tensile_stress/MPa)
    set_cell_value(book, f"{sheet_name}!bolt.E", fseg.bolt.elastic_modulus/GPa)
    set_cell_value(book, f"{sheet_name}!bolt.pretension", fseg.Fv/kN)
    set_cell_value(book, f"{sheet_name}!Do", fseg.Do/mm)
    set_cell_value(book, f"{sheet_name}!Dw", fseg.Dw/mm)
    set_cell_value(book, f"{sheet_name}!gap.h", fseg.gap_height/mm)
    set_cell_value(book, f"{sheet_name}!gap.L", fseg.gap_angle * fseg.R/mm)
    set_cell_value(book, f"{sheet_name}!E_mod", fseg.E/GPa)
    set_cell_value(book, f"{sheet_name}!G_mod", fseg.G/GPa)

    # Polynomial Data Series
    Z = np.linspace(1.05*fseg._compressive_force_polynomial.domain[0], 1.05*fseg.shell_force_at_tensile_ULS, 100)
    Fs = fseg.bolt_axial_force(Z)
    Ms = fseg.bolt_bending_moment(Z)
    set_cell_value(book, f"{sheet_name}!dataseries.Z", Z/kN)
    set_cell_value(book, f"{sheet_name}!dataseries.Fs", Fs/kN)
    set_cell_value(book, f"{sheet_name}!dataseries.Ms", Ms)

    # Polynomial Data Points
    Z1 = fseg.shell_force_at_rest
    Fs1 = fseg.bolt_force_at_rest
    Ms1 = fseg.bolt_moment_at_rest

    Z2 = fseg.shell_force_at_tensile_ULS
    Fs2 = fseg.bolt_force_at_tensile_ULS
    Ms2 = fseg.bolt_moment_at_tensile_ULS

    Z3 = fseg.shell_force_at_small_displacement
    Fs3 = fseg.bolt_force_at_small_displacement
    Ms3 = fseg.bolt_moment_at_small_displacement

    Z4 = fseg._compressive_force_polynomial.domain[0]
    Fs4 = fseg.bolt_axial_force(Z4)
    Ms4 = fseg.bolt_bending_moment(Z4)

    set_cell_value(book, f"{sheet_name}!point1.Z", Z1/kN)
    set_cell_value(book, f"{sheet_name}!point2.Z", Z2/kN)
    set_cell_value(book, f"{sheet_name}!point3.Z", Z3/kN)
    set_cell_value(book, f"{sheet_name}!point4.Z", Z4/kN)

    set_cell_value(book, f"{sheet_name}!point1.Fs", Fs1/kN)
    set_cell_value(book, f"{sheet_name}!point2.Fs", Fs2/kN)
    set_cell_value(book, f"{sheet_name}!point3.Fs", Fs3/kN)
    set_cell_value(book, f"{sheet_name}!point4.Fs", Fs4/kN)

    set_cell_value(book, f"{sheet_name}!point1.Ms", Ms1)
    set_cell_value(book, f"{sheet_name}!point2.Ms", Ms2)
    set_cell_value(book, f"{sheet_name}!point3.Ms", Ms3)
    set_cell_value(book, f"{sheet_name}!point4.Ms", Ms4)

    # Polynomial Data Points
    set_cell_value(book, f"{sheet_name}!Fs_coeff.tens", fseg._tensile_force_polynomial.coef)
    set_cell_value(book, f"{sheet_name}!Fs_coeff.comp", fseg._compressive_force_polynomial.coef)
    set_cell_value(book, f"{sheet_name}!Fs_coeff.const", [Fs4, 0, 0])

    set_cell_value(book, f"{sheet_name}!Ms_coeff.tens", fseg._tensile_moment_polynomial.coef)
    set_cell_value(book, f"{sheet_name}!Ms_coeff.comp", fseg._compressive_moment_polynomial.coef)
    set_cell_value(book, f"{sheet_name}!Ms_coeff.const", [Ms4, 0, 0])

    # Miscellaneous Model Data
    set_cell_value(book, f"{sheet_name}!gap.angle", fseg.gap_angle/deg)
    set_cell_value(book, f"{sheet_name}!bolt.As", fseg.bolt.tensile_cross_section_area/(mm**2))
    set_cell_value(book, f"{sheet_name}!a_prime", (fseg.b / (fseg._prying_lever_ratio - 1))/mm)
    set_cell_value(book, f"{sheet_name}!Z0", fseg._ideal_shell_force_at_tensile_ULS/kN)
    set_cell_value(book, f"{sheet_name}!Z2_tilde", fseg._cantilever_shell_force_at_tensile_ULS/kN)
    set_cell_value(book, f"{sheet_name}!bolt.axial_stiffness", fseg._bolt_axial_stiffness/(kN/mm))
    set_cell_value(book, f"{sheet_name}!bolt.bending_stiffness", fseg._bolt_bending_stiffness/kN)
    set_cell_value(book, f"{sheet_name}!clamped_parts_stiffness", fseg._flange_axial_stiffness/(kN/mm))
    set_cell_value(book, f"{sheet_name}!gap.stiffness", fseg._gap_stiffness/1e6)
    set_cell_value(book, f"{sheet_name}!DZ_gap", fseg.shell_force_at_closed_gap/kN)
    set_cell_value(book, f"{sheet_name}!stiffness_correction_factor", fseg._stiffness_correction_factor)
    set_cell_value(book, f"{sheet_name}!polynomial_initial_slope", fseg._polynomial_initial_slope)
    set_cell_value(book, f"{sheet_name}!true_force_initial_slope", fseg._tensile_force_polynomial.deriv()(Z1))
    set_cell_value(book, f"{sheet_name}!true_moment_initial_slope", fseg._tensile_moment_polynomial.deriv()(Z1)*1000)


wb = Book(os.path.join(os.path.dirname(__file__), "BnB_3p228mm-Results.xlsx"))

print("")
print("Flange Segment Model with 30 deg gap width")
print("------------------------------------------")
fseg_30deg  = create_flange_segment( 30*deg)
flange_segment_model_to_excel(wb, "Gap30deg", fseg_30deg)

print("")
print("Flange Segment Model with 60 deg gap width")
print("------------------------------------------")
fseg_60deg  = create_flange_segment( 60*deg)
flange_segment_model_to_excel(wb, "Gap60deg", fseg_60deg)

print("")
print("Flange Segment Model with 90 deg gap width")
print("------------------------------------------")
fseg_90deg  = create_flange_segment( 90*deg)
flange_segment_model_to_excel(wb, "Gap90deg", fseg_90deg)

print("")
print("Flange Segment Model with 120 deg gap width")
print("-------------------------------------------")
fseg_120deg = create_flange_segment(120*deg)
flange_segment_model_to_excel(wb, "Gap120deg", fseg_120deg)
