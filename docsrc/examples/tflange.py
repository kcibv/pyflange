# PyFlange tutorial: deterministic fatigue analysis of a T-flange.
#
# This script is explained step by step in the tutorial
# https://kcibv.github.io/pyflange/latest/tutorials/tflange/
#
# All values are in SI units: m, N, Pa, rad.

from math import pi
import numpy as np

mm  = 0.001     # m
kN  = 1000      # N
MPa = 1e6       # Pa
deg = pi/180    # rad


# --8<-- [start:segment]
from pyflange.bolts import StandardMetricBolt, ISOHexNut, standard_hole_diameter
from pyflange.flangesegments import PolynomialTFlangeSegment, Gap, shell_stiffness
from pyflange.stats import gap_height_distribution

D = 7500*mm                 # tower outer diameter
n_pairs = 200               # number of bolt pairs (one inside and one outside the shell)
s = 90*mm                   # shell thickness

gap_angle = 30*deg
gap_dist = gap_height_distribution(D, 0.0014, gap_angle * D/2)

fseg = PolynomialTFlangeSegment(
    a = 62.5*mm,            # flange edge to bolt axis
    b = 111*mm,             # bolt axis to shell center line
    s = s,
    t = 120*mm,             # flange thickness
    R = D/2,                # shell outer radius
    central_angle = 2*pi / n_pairs,

    Zg = -105.8*kN,         # dead weight per bolt pair (compression)

    bolt = StandardMetricBolt("M48", "10.9", shank_length=150*mm,
                              shank_diameter_ratio="pitch", stud=True),
    Fv = 928*kN,            # preload of each bolt, before losses
    Do = standard_hole_diameter("M48"),
    washer = None,
    nut = ISOHexNut("M48"),

    gap = Gap(height = gap_dist.ppf(0.95), angle = gap_angle),
    tilt_angle = 0.0,       # required: no taper effect for this flange

    k_shell = shell_stiffness(D/2, s, gap_angle))
# --8<-- [end:segment]


# --8<-- [start:validate]
mode, Zu = fseg.failure_mode(fy_sh=355*MPa, fy_fl=355*MPa)
print(f"Governing failure mode: {mode}")
print(f"Ultimate shell forces A, B, D, E [kN]: {[round(Z/kN) for Z in Zu]}")

fseg.validate(fy_sh=355*MPa, fy_fl=355*MPa)    # T-flanges: failure mode A must govern
# --8<-- [end:validate]


# --8<-- [start:transfer]
Z = np.array([0, 250, 500, 1000]) * kN      # shell force per bolt pair
for Zi, Fs, Ms in zip(Z, fseg.bolt_axial_force(Z), fseg.bolt_bending_moment(Z)):
    print(f"Z = {Zi/kN:6.0f} kN   Fs = {Fs/kN:6.1f} kN   Ms = {Ms:6.0f} Nm")
# --8<-- [end:transfer]


# --8<-- [start:fatigue]
from pyflange.fatigue import MarkovMatrix, BoltFatigueAnalysis

mkvm = MarkovMatrix(
    range  = np.array([ 20,  40,  80, 120, 160, 200, 240, 300, 400, 500]) * 1e6,   # N·m
    mean   = np.array([ 80,  80,  80,  80,  80,  80,  80,  80,  80,  80]) * 1e6,   # N·m
    cycles = np.array([3e8, 1e8, 3e7, 8e6, 2e6, 5e5, 1.5e5, 4e4, 8e3, 1e3]),
    duration = 25)                                                               # years

fatigue = BoltFatigueAnalysis(fseg, mkvm)
print(f"Damage:       {fatigue.damage:.3f}")
print(f"Fatigue life: {fatigue.fatigue_life:.0f} years")
# --8<-- [end:fatigue]
