import sys
from pyflange.utils import Logger, log_data
logger = Logger(__name__)

from pyflange.flangesegments import PolynomialLFlangeSegment, shell_stiffness, Gap
from pyflange.bolts import MetricBolt, HexNut, FlatWasher
from pyflange.gap import gap_height_distribution

from math import *
import numpy as np

from workbook import *



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


# Bolts

M48 = MetricBolt(
    nominal_diameter = 48*mm,
    thread_pitch = 5*mm,
    shank_length = 238*mm,
    yield_stress = 900*MPa,
    ultimate_tensile_stress = 1040*MPa,
    stud = False)

M48_hex_nut = HexNut(
    nominal_diameter = 48*mm,
    thickness = 38*mm,
    inscribed_diameter = 75*mm,
    circumscribed_diameter = 82.6*mm,
    bearing_diameter = 69.4*mm
)

M80 = MetricBolt(
    nominal_diameter = 80*mm,
    thread_pitch = 6*mm,
    shank_diameter_ratio = 76.1/80,
    shank_length = 270*mm,
    yield_stress = 900*MPa,
    ultimate_tensile_stress = 1000*MPa,
    stud = True)

M80_hex_nut = HexNut(
    nominal_diameter = 80*mm,
    thickness = 64*mm,
    inscribed_diameter = 115*mm,
    circumscribed_diameter = 127.5*mm,
    bearing_diameter = 140*mm
)


# Polinomial Segment Models

def create_D4600_flange_segment (gap_angle, gap_shape_factor=1.0, tilt_angle=0, interp_shell_stiff=False):

    D = 4600*mm
    t_sh = 30*mm
    n = 144 # number of bolts
    gap_length = gap_angle * D/2
    gap = gap_height_distribution(D, 0.0014, gap_length)

    k_mean = gap.mean()
    COV_k = gap.std() / k_mean

    fseg = PolynomialLFlangeSegment(

        a = 160*mm,           # distance between inner face of the flange and center of the bolt hole
        b = 80*mm,           # distance between center of the bolt hole and center-line of the shell
        s = t_sh,               # shell thickness
        t = 125*mm,           # flange thickness
        R = D/2,                # shell outer curvature radius
        central_angle = 2*pi/n, # angle subtented by the flange segment arc

        Zg = -30.139*kN,     # load applied to the flange segment shell at rest
                                # (normally dead weight of tower + RNA, divided by the number of bolts)

        bolt = M48,
        Fv = 910*kN,       # applied bolt preload

        Do = 52*mm,         # bolt hole diameter
        washer = FlatWasher(outer_diameter=92*mm, inner_diameter=49.4*mm, thickness=8*mm),   # washer
        nut = M48_hex_nut,

        gap = Gap(height = gap.ppf(0.95),             # maximum longitudinal gap height
                  angle = gap_angle,                  # longitudinal gap length
                  shape_factor = gap_shape_factor),   # scaling factor accounting for the gap shape

        tilt_angle = tilt_angle,    # flange tilt angle

        s_ratio = 30/30,    # ratio of bottom shell thickness over tower shell thickness

        k_shell = shell_stiffness(D/2, (t_sh + 30*mm)/2, gap_angle) if interp_shell_stiff else None
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    log_data(fseg, k_mean=k_mean, COV_k=COV_k)

    return fseg


def create_D7500_flange_segment (gap_angle, gap_shape_factor=1.0, tilt_angle=0, interp_shell_stiff=False):

    D = 7500*mm
    t_sh = 72*mm
    n = 120 # number of bolts
    gap_length = gap_angle * D/2
    gap = gap_height_distribution(D, 0.0014, gap_length)

    k_mean = gap.mean()
    COV_k = gap.std() / k_mean

    fseg = PolynomialLFlangeSegment(

        a = 232.5*mm,           # distance between inner face of the flange and center of the bolt hole
        b = 166.5*mm,           # distance between center of the bolt hole and center-line of the shell
        s = t_sh,               # shell thickness
        t = 200.0*mm,           # flange thickness
        R = D/2,                # shell outer curvature radius
        central_angle = 2*pi/n, # angle subtented by the flange segment arc

        Zg = -14795*kN / n,     # load applied to the flange segment shell at rest
                                # (normally dead weight of tower + RNA, divided by the number of bolts)

        bolt = M80,
        Fv = 2800*kN,       # applied bolt preload

        Do = 86*mm,         # bolt hole diameter
        washer = None,      # no washer
        nut = M80_hex_nut,  # bolt nut

        gap = Gap(height = gap.ppf(0.95),             # maximum longitudinal gap height
                  angle = gap_angle,                  # longitudinal gap length
                  shape_factor = gap_shape_factor),    # scaling factor accounting for the gap shape

        tilt_angle = tilt_angle,    # flange tilt angle

        s_ratio = 100/72,    # ratio of bottom shell thickness over tower shell thickness

        k_shell = shell_stiffness(D/2, (t_sh + 100*mm)/2, gap_angle) if interp_shell_stiff else None
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    log_data(fseg, k_mean=k_mean, COV_k=COV_k)

    return fseg



def damage (fseg, flange_markov_matrix):
    from pyflange.fatigue import BoltFatigueCurve
    from pyflange.flangesegments import bolt_markov_matrix
    kb = max(0.5, 0.5 + 0.5*log(fseg.bolt.nominal_diameter/0.036) / log(150/36))
    snc = BoltFatigueCurve(fseg.bolt.nominal_diameter)
    bolt_mm = bolt_markov_matrix(fseg, flange_markov_matrix, bending_factor=kb)
    return snc.cumulated_damage(bolt_mm)


import sys
params = sys.argv[1:]


if len(params) == 0 or "-h" in params:
    import os
    print("Generate the PyFlange validation report cases." )
    print("")
    print("Usage:")
    print("    gen-cases.py <case-list> [-h | -i]")
    print("")
    print("Parameters:")
    print("    <case-list>  is the list of the case numbers to be generated. One or more of:")
    print("                 1 (D = 7500 mm with sinusoidal gap shape and no flange tilt.)")
    print("                 2 (D = 7500 mm with sinusoidal gap shape and 1 deg flange tilt")
    print("                 3 (D = 7500 mm gap shape factor 1.2 and no flange tilt")
    print("                 4 (D = 4600 mm with sinusoidal gap shape and no flange tilt")
    print("")
    print("    -i           uses the interpolated shell stiffness instead of the simplified formula")
    print("    -h           prints this help message")
    print("")
    print("Example:")
    print("    gen-cases.py 1 3 4 -i")
    print("")



if "1" in params:
    print("\nEvaluating D7500 Flange Segment Model with sinusoidal gap shape and no flange tilt ...")
    wb1 = open_workbook(f"Case_1_D7500_L_Tilt-0p00deg_ShapeFactor-1p00_Fv2800kN_ShellStiff-{'Interp' if '-i' in params else 'Simpl'}.xlsx")
    flange_markov_matrix = load_markov_matrix(wb1, "SGRE!markov_matrix")

    print("... with 30 deg gap width")
    fseg_30deg  = create_D7500_flange_segment( 30*deg, interp_shell_stiff = "-i" in params)
    fseg_30deg._damage = damage(fseg_30deg, flange_markov_matrix)
    flangesegment_to_excel(wb1, "PyFlange_Gap30deg", fseg_30deg)

    print("... with 60 deg gap width")
    fseg_60deg  = create_D7500_flange_segment( 60*deg, interp_shell_stiff = "-i" in params)
    fseg_60deg._damage = damage(fseg_60deg, flange_markov_matrix)
    flangesegment_to_excel(wb1, "PyFlange_Gap60deg", fseg_60deg)

    print("... with 90 deg gap width")
    fseg_90deg  = create_D7500_flange_segment( 90*deg, interp_shell_stiff = "-i" in params)
    fseg_90deg._damage = damage(fseg_90deg, flange_markov_matrix)
    flangesegment_to_excel(wb1, "PyFlange_Gap90deg", fseg_90deg)

    print("... with 120 deg gap width")
    fseg_120deg = create_D7500_flange_segment(120*deg, interp_shell_stiff = "-i" in params)
    fseg_120deg._damage = damage(fseg_120deg, flange_markov_matrix)
    flangesegment_to_excel(wb1, "PyFlange_Gap120deg", fseg_120deg)



if "2" in params:
    print("\nEvaluating D7500 Flange Segment Model with sinusoidal gap shape and 1 deg flange tilt ...")
    wb2 = open_workbook(f"Case_2_D7500_L_Tilt-0p25deg_ShapeFactor-1p00_Fv2800kN_ShellStiff-{'Interp' if '-i' in params else 'Simpl'}.xlsx")
    flange_markov_matrix = load_markov_matrix(wb2, "SGRE!markov_matrix")

    print("... with 30 deg gap width")
    fseg_30deg  = create_D7500_flange_segment( 30*deg, tilt_angle=0.25*deg, interp_shell_stiff = "-i" in params)
    fseg_30deg._damage = damage(fseg_30deg, flange_markov_matrix)
    flangesegment_to_excel(wb2, "PyFlange_Gap30deg", fseg_30deg)

    print("... with 60 deg gap width")
    fseg_60deg  = create_D7500_flange_segment( 60*deg, tilt_angle=0.25*deg, interp_shell_stiff = "-i" in params)
    fseg_60deg._damage = damage(fseg_60deg, flange_markov_matrix)
    flangesegment_to_excel(wb2, "PyFlange_Gap60deg", fseg_60deg)

    print("... with 90 deg gap width")
    fseg_90deg  = create_D7500_flange_segment( 90*deg, tilt_angle=0.25*deg, interp_shell_stiff = "-i" in params)
    fseg_90deg._damage = damage(fseg_90deg, flange_markov_matrix)
    flangesegment_to_excel(wb2, "PyFlange_Gap90deg", fseg_90deg)

    print("... with 120 deg gap width")
    fseg_120deg = create_D7500_flange_segment(120*deg, tilt_angle=0.25*deg, interp_shell_stiff = "-i" in params)
    fseg_120deg._damage = damage(fseg_120deg, flange_markov_matrix)
    flangesegment_to_excel(wb2, "PyFlange_Gap120deg", fseg_120deg)



if "3" in params:
    print("\nEvaluating D7500 Flange Segment Model with gap shape factor 1.2 and no flange tilt ...")
    wb3 = open_workbook(f"Case_3_D7500_L_Tilt-0p00deg_ShapeFactor-1p20_Fv2800kN_ShellStiff-{'Interp' if '-i' in params else 'Simpl'}.xlsx")
    flange_markov_matrix = load_markov_matrix(wb3, "SGRE!markov_matrix")

    print("... with 30 deg gap width")
    fseg_30deg  = create_D7500_flange_segment( 30*deg, 1.2, interp_shell_stiff = "-i" in params)
    fseg_30deg._damage = damage(fseg_30deg, flange_markov_matrix)
    flangesegment_to_excel(wb3, "PyFlange_Gap30deg", fseg_30deg)

    print("... with 60 deg gap width")
    fseg_60deg  = create_D7500_flange_segment( 60*deg, 1.2, interp_shell_stiff = "-i" in params)
    fseg_60deg._damage = damage(fseg_60deg, flange_markov_matrix)
    flangesegment_to_excel(wb3, "PyFlange_Gap60deg", fseg_60deg)

    print("... with 90 deg gap width")
    fseg_90deg  = create_D7500_flange_segment( 90*deg, 1.2, interp_shell_stiff = "-i" in params)
    fseg_90deg._damage = damage(fseg_90deg, flange_markov_matrix)
    flangesegment_to_excel(wb3, "PyFlange_Gap90deg", fseg_90deg)

    print("... with 120 deg gap width")
    fseg_120deg = create_D7500_flange_segment(120*deg, 1.2, interp_shell_stiff = "-i" in params)
    fseg_120deg._damage = damage(fseg_120deg, flange_markov_matrix)
    flangesegment_to_excel(wb3, "PyFlange_Gap120deg", fseg_120deg)



if "4" in params:
    print("\nEvaluating D4600 Flange Segment Model with sinusoidal gap shape and no flange tilt ...")
    wb4 = open_workbook(f"Case_4_D4600_L_Tilt-0p00deg_ShapeFactor-1p00_Fv910kN_ShellStiff-{'Interp' if '-i' in params else 'Simpl'}.xlsx")
    flange_markov_matrix = load_markov_matrix(wb4, "SGRE!markov_matrix")

    print("... with 30 deg gap width")
    fseg_30deg  = create_D4600_flange_segment( 30*deg, interp_shell_stiff = "-i" in params)
    fseg_30deg._damage = damage(fseg_30deg, flange_markov_matrix)
    flangesegment_to_excel(wb4, "PyFlange_Gap30deg", fseg_30deg)

    print("... with 60 deg gap width")
    fseg_60deg  = create_D4600_flange_segment( 60*deg, interp_shell_stiff = "-i" in params)
    fseg_60deg._damage = damage(fseg_60deg, flange_markov_matrix)
    flangesegment_to_excel(wb4, "PyFlange_Gap60deg", fseg_60deg)

    print("... with 90 deg gap width")
    fseg_90deg  = create_D4600_flange_segment( 90*deg, interp_shell_stiff = "-i" in params)
    fseg_90deg._damage = damage(fseg_90deg, flange_markov_matrix)
    flangesegment_to_excel(wb4, "PyFlange_Gap90deg", fseg_90deg)

    print("... with 120 deg gap width")
    fseg_120deg = create_D4600_flange_segment(120*deg, interp_shell_stiff = "-i" in params)
    fseg_120deg._damage = damage(fseg_120deg, flange_markov_matrix)
    flangesegment_to_excel(wb4, "PyFlange_Gap120deg", fseg_120deg)
