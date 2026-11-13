from math import pi
from pyflange.bolts import StandardMetricBolt, RoundNut
from pyflange.flangesegments import PolynomialLFlangeSegment, Gap

# Units: meters, newtons, radians
mm = 0.001
kN = 1000
deg = pi/180

# The flange segment: one bolt of an L-flange with 156 M80 bolts
fseg = PolynomialLFlangeSegment(
    a = 150*mm,                     # inner face of the flange to bolt axis
    b = 122*mm,                     # bolt axis to shell center line
    s = 54*mm,                      # shell thickness
    t = 172*mm,                     # flange thickness
    R = 4000*mm,                    # shell outer radius
    central_angle = 2*pi/156,       # 156 bolts
    Zg = -115*kN,                   # dead weight per bolt
    bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm, stud=True),
    Fv = 2876*kN,                   # bolt preload
    Do = 86*mm,                     # bolt hole diameter
    washer = None,
    nut = RoundNut("M80"),
    gap = Gap(height=0.89*mm, angle=30*deg))

# Print the bolt force for some values of the shell force
for Z in [-200*kN, 0*kN, 200*kN, 400*kN, 600*kN, 800*kN, 1000*kN]:
    Fs = fseg.bolt_axial_force(Z)
    print(f"Shell force {Z/kN:6.0f} kN  ->  bolt force {Fs/kN:7.1f} kN")
