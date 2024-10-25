from pyflange.logger import Logger, log_data
logger = Logger(__name__)

from pyflange.flangesegments import PolynomialLFlangeSegment, shell_stiffness
from pyflange.bolts import MetricBolt, HexNut
from pyflange.gap import gap_height_distribution

from math import pi
import numpy as np

from workbook import open_workbook, flangesegment_to_excel



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


M80_hex_nut = HexNut(
    nominal_diameter = 80*mm,
    thickness = 64*mm,
    inscribed_diameter = 115*mm,
    circumscribed_diameter = 127.5*mm,
    bearing_diameter = 140*mm
)


# Polinomial Segment Model

def create_flange_segment (gap_angle, gap_shape_factor=1.0, tilt_angle=0):

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

        gap_height = gap.ppf(0.95),             # maximum longitudinal gap height
        gap_angle = gap_angle,                  # longitudinal gap length
        gap_shape_factor = gap_shape_factor,    # scaling factor accounting for the gap shape

        tilt_angle = tilt_angle,    # flange tilt angle

        s_ratio = 100/72    # ratio of bottom shell thickness over tower shell thickness

        #k_shell = shell_stiffness(D/2, (t_sh + 100*mm)/2, gap_angle)
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    log_data(fseg, k_mean=k_mean, COV_k=COV_k)

    return fseg



print("\nEvaluating Flange Segment Model with sinusoidal gap shape and no flange tilt ...")
wb = open_workbook("Case-L1_D7500_Tilt-0p00deg_ShapeFactor-1p00_Fv2800kN.xlsx")

print("... with 30 deg gap width")
fseg_30deg  = create_flange_segment( 30*deg)
flangesegment_to_excel(wb, "PyFlange_Gap30deg", fseg_30deg)

print("... with 60 deg gap width")
fseg_60deg  = create_flange_segment( 60*deg)
flangesegment_to_excel(wb, "PyFlange_Gap60deg", fseg_60deg)

print("... with 90 deg gap width")
fseg_90deg  = create_flange_segment( 90*deg)
flangesegment_to_excel(wb, "PyFlange_Gap90deg", fseg_90deg)

print("... with 120 deg gap width")
fseg_120deg = create_flange_segment(120*deg)
flangesegment_to_excel(wb, "PyFlange_Gap120deg", fseg_120deg)




print("\nEvaluating Flange Segment Model with sinusoidal gap shape and 1 deg flange tilt ...")
wb_tt = open_workbook("Case-L2_D7500_Tilt-0p25deg_ShapeFactor-1p00_Fv2800kN.xlsx")

print("... with 30 deg gap width")
fseg_30deg_tt  = create_flange_segment( 30*deg, tilt_angle=0.25*deg)
flangesegment_to_excel(wb_tt, "PyFlange_Gap30deg", fseg_30deg_tt)

print("... with 60 deg gap width")
fseg_60deg_tt  = create_flange_segment( 60*deg, tilt_angle=0.25*deg)
flangesegment_to_excel(wb_tt, "PyFlange_Gap60deg", fseg_60deg_tt)

print("... with 90 deg gap width")
fseg_90deg_tt  = create_flange_segment( 90*deg, tilt_angle=0.25*deg)
flangesegment_to_excel(wb_tt, "PyFlange_Gap90deg", fseg_90deg_tt)

print("... with 120 deg gap width")
fseg_120deg_tt = create_flange_segment(120*deg, tilt_angle=0.25*deg)
flangesegment_to_excel(wb_tt, "PyFlange_Gap120deg", fseg_120deg_tt)



print("\nEvaluating Flange Segment Model with gap shape factor 1.2 and no flange tilt ...")
wb_sf = open_workbook("Case-L3_D7500_Tilt-0p00deg_ShapeFactor-1p20_Fv2800kN.xlsx")

print("... with 30 deg gap width")
fseg_30deg_sf  = create_flange_segment( 30*deg, 1.2)
flangesegment_to_excel(wb_sf, "PyFlange_Gap30deg", fseg_30deg_sf)

print("... with 60 deg gap width")
fseg_60deg_sf  = create_flange_segment( 60*deg, 1.2)
flangesegment_to_excel(wb_sf, "PyFlange_Gap60deg", fseg_60deg_sf)

print("... with 90 deg gap width")
fseg_90deg_sf  = create_flange_segment( 90*deg, 1.2)
flangesegment_to_excel(wb_sf, "PyFlange_Gap90deg", fseg_90deg_sf)

print("... with 120 deg gap width")
fseg_120deg_sf = create_flange_segment(120*deg, 1.2)
flangesegment_to_excel(wb_sf, "PyFlange_Gap120deg", fseg_120deg_sf)
