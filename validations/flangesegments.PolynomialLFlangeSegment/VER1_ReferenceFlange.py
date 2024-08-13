from pyflange.logger import Logger, log_data
logger = Logger(__name__)

from pyflange.flangesegments import PolynomialLFlangeSegment
from pyflange.bolts import MetricBolt
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

M48 = MetricBolt(
    nominal_diameter = 48*mm,
    thread_pitch = 5*mm,
    shank_length = 238*mm,
    yield_stress = 900*MPa,
    ultimate_tensile_stress = 1040*MPa,
    stud = True)





# Polinomial Segment Model

def create_flange_segment (gap_angle, gap_shape_factor=1.0, tilt_angle=0):

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
        Fv = 809*kN,       # applied bolt preload

        Do = 52*mm,         # bolt hole diameter
        Dw = 92*mm,        # washer diameter

        gap_height = gap.ppf(0.95),             # maximum longitudinal gap height
        gap_angle = gap_angle,                  # longitudinal gap length
        gap_shape_factor = gap_shape_factor,    # scaling factor accounting for the gap shape

        tilt_angle = tilt_angle,    # flange tilt angle

        s_ratio = 30/30    # ratio of bottom shell thickness over tower shell thickness
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    log_data(fseg, k_mean=k_mean, COV_k=COV_k)

    return fseg







print("\nEvaluating Flange Segment Model with sinusoidal gap shape and no flange tilt ...")
wb = open_workbook("VER1_ReferenceFlange-Results.xlsx")

print("... with 30 deg gap width")
fseg_30deg  = create_flange_segment( 30*deg)
flangesegment_to_excel(wb, "Gap30deg", fseg_30deg)

print("... with 60 deg gap width")
fseg_60deg  = create_flange_segment( 60*deg)
flangesegment_to_excel(wb, "Gap60deg", fseg_60deg)

print("... with 90 deg gap width")
fseg_90deg  = create_flange_segment( 90*deg)
flangesegment_to_excel(wb, "Gap90deg", fseg_90deg)

print("... with 120 deg gap width")
fseg_120deg = create_flange_segment(120*deg)
flangesegment_to_excel(wb, "Gap120deg", fseg_120deg)
