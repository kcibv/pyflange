# PyFlange tutorial: design gap angles and comparison of AMD1 and ED2.
#
# This script is explained step by step in the tutorial
# https://kcibv.github.io/pyflange/latest/tutorials/design-gaps/
#
# All values are in SI units: m, N, Pa, rad.

from math import pi
import numpy as np

mm  = 0.001     # m
kN  = 1000      # N
deg = pi/180    # rad


# --8<-- [start:common]
from pyflange.bolts import StandardMetricBolt, RoundNut, standard_hole_diameter
from pyflange.fatigue import MarkovMatrix, BoltFatigueAnalysis

D = 8000*mm                 # tower outer diameter
R = D/2
s = 54*mm
u_tol = 0.0014              # flatness tolerance, 1.4 mm/m
n_bolts = 156

# Parameters shared by all the flange segments of this study
common = dict(
    a = 150*mm, b = 122*mm, s = s, t = 172*mm, R = R,
    central_angle = 2*pi / n_bolts,
    Zg = -18044*kN / n_bolts,
    bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm,
                              shank_diameter_ratio="pitch", stud=True),
    Fv = 2876*kN,
    Do = standard_hole_diameter("M80"),
    washer = None,
    nut = RoundNut("M80"),
    tilt_angle = 0.1*deg,
    PPL = 0.05, DFT = 300e-6, settlement_factor = 0.5)

mkvm = MarkovMatrix(
    range  = np.array([ 25,  50, 100, 150, 200, 250,  300, 375, 500, 625]) * 1e6,
    mean   = np.array([100, 100, 100, 100, 100, 100,  100, 100, 100, 100]) * 1e6,
    cycles = np.array([3e8, 1e8, 3e7, 8e6, 2e6, 5e5, 1.5e5, 4e4, 8e3, 1e3]),
    duration = 30)
# --8<-- [end:common]


# --8<-- [start:amd1]
from pyflange.flangesegments import PolynomialLFlangeSegment, Gap, shell_stiffness
from pyflange.stats import gap_height_distribution

print("AMD1: design gap = 95% quantile of the log-normal distribution")
amd1_damage = {}
for angle_deg in (30, 60, 90, 120):
    angle = angle_deg * deg
    k_design = gap_height_distribution(D, u_tol, angle * D/2).ppf(0.95)
    fseg = PolynomialLFlangeSegment(**common,
                                    gap = Gap(height=k_design, angle=angle),
                                    k_shell = shell_stiffness(R, s, angle))
    amd1_damage[angle_deg] = BoltFatigueAnalysis(fseg, mkvm).damage
    print(f"  gap angle {angle_deg:3d} deg   gap height {k_design/mm:4.2f} mm   damage {amd1_damage[angle_deg]:.3f}")

governing = max(amd1_damage, key=amd1_damage.get)
print(f"  governing gap angle: {governing} deg")
# --8<-- [end:amd1]


# --8<-- [start:ed2]
from pyflange.flangesegments import PolynomialLFlangeSegmentED2
from pyflange.stats import gap_height_distribution_ED2

print("ED2: design gap = 75% quantile of the Weibull distribution")
for angle_deg in (30, 60, 90, 120):
    angle = angle_deg * deg
    k_design = gap_height_distribution_ED2(D, u_tol, angle * D/2).ppf(0.75)
    fseg = PolynomialLFlangeSegmentED2(**common,
                                       gap = Gap(height=k_design, angle=angle),
                                       k_shell = shell_stiffness(R, s, angle))
    damage = BoltFatigueAnalysis(fseg, mkvm).damage
    print(f"  gap angle {angle_deg:3d} deg   gap height {k_design/mm:4.2f} mm   damage {damage:.3f}")
# --8<-- [end:ed2]
