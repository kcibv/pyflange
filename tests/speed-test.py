from datetime import datetime

from pyflange.flangesegments import PolynomialLFlangeSegment
from pyflange.bolts import StandardMetricBolt, ISOFlatWasher, ISOHexNut
from pyflange.gap import gap_height_distribution

from math import pi
import numpy as np


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


def log (message):
    ''' Logs the given message, preceded by a timestamp '''
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
    print("[{}] {}".format(timestamp, message))


M80 = StandardMetricBolt("M80", "10.9", 
    shank_diameter_ratio = 76.1/80,
    shank_length = 270*mm,
    stud = True)


def create_flange_segment (gap_angle):

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
        s = t_sh,           # shell thickness
        t = 200.0*mm,           # flange thickness
        central_angle = 2*pi/n,   # shell arc angle
        R = D/2,          # shell outer curvature radius

        Zg = -14795*kN / n, # load applied to the flange segment shell at rest
                                                # (normally dead weight of tower + RNA, divided by the number of bolts)

        bolt = M80,
        Fv = 2876*kN,                            # applied bolt preload

        Do = 86*mm,     # bolt hole diameter
        washer = ISOFlatWasher("M80"),    # washer
        nut = ISOHexNut("M80"), 

        gap_height = gap.ppf(0.95),   # maximum longitudinal gap height
        gap_angle = gap_angle,  # longitudinal gap length

        s_ratio = 102/72        # ratio of bottom shell thickness over tower shell thickness
    )

    # Assert that failure mode is B.
    #fseg.validate(335*MPa, 285*MPa)

    return fseg


def test_speed ():
    n = 10000
    Z = np.linspace(-1000, 1000, 1000)
    log(f"Starting PyFlange speed test for {n} iterations...")
    for i in range(n):
        fseg = create_flange_segment(pi/2)
        Fs = fseg.bolt_axial_force(Z)
        Ms = fseg.bolt_bending_moment(Z)
        if i % 1000 == 0 and i != 0:
            log(f"... performed {i} iterations")
    log(f"... DONE")

if __name__ == "__main__":
    test_speed()
