# PyFlange tutorial: deterministic fatigue analysis of an L-flange.
#
# This script is explained step by step in the tutorial
# https://kcibv.github.io/pyflange/latest/tutorials/lflange/
#
# All values are in SI units: m, N, Pa, rad.

# --8<-- [start:units]
from math import pi

mm  = 0.001     # m
kN  = 1000      # N
MPa = 1e6       # Pa
deg = pi/180    # rad
# --8<-- [end:units]


# --8<-- [start:bolt]
from pyflange.bolts import StandardMetricBolt, RoundNut, standard_hole_diameter

bolt = StandardMetricBolt("M80", "10.9",
                          shank_length = 160*mm,
                          shank_diameter_ratio = "pitch",   # reduced shank
                          stud = True)
nut = RoundNut("M80")
Do = standard_hole_diameter("M80")      # 86 mm bolt hole
# --8<-- [end:bolt]


# --8<-- [start:gap]
from pyflange.stats import gap_height_distribution
from pyflange.flangesegments import Gap, shell_stiffness

D = 8000*mm                 # flange (tower) outer diameter
u_tol = 1.4*mm / 1.0        # flatness tolerance: 1.4 mm/m
gap_angle = 30*deg          # gap angle
gap_length = gap_angle * D/2

gap_dist = gap_height_distribution(D, u_tol, gap_length)
gap = Gap(height = gap_dist.ppf(0.95),   # design gap: 95% quantile
          angle = gap_angle)
# --8<-- [end:gap]


# --8<-- [start:segment]
from pyflange.flangesegments import PolynomialLFlangeSegment

n_bolts = 156               # number of bolts in the flange
R = D/2                     # shell outer radius
s = 54*mm                   # shell thickness

fseg = PolynomialLFlangeSegment(
    a = 150*mm,             # inner face of the flange to bolt axis
    b = 122*mm,             # bolt axis to shell center line
    s = s,
    t = 172*mm,             # flange thickness
    R = R,
    central_angle = 2*pi / n_bolts,

    Zg = -18044*kN / n_bolts,           # dead weight per bolt (compression)

    bolt = bolt,
    Fv = 2876*kN,                       # bolt preload before losses
    Do = Do,
    washer = None,
    nut = nut,

    gap = gap,
    tilt_angle = 0.1*deg,               # inside flange inclination

    k_shell = shell_stiffness(R, s, gap_angle),   # shell stiffness from FEA data
    PPL = 0.05,                         # 5% preload loss due to plastic strain
    DFT = 300e-6,                       # 300 um coating thickness
    settlement_factor = 0.5)            # tension tightening with re-tightening
# --8<-- [end:segment]


# --8<-- [start:validate]
mode, Zu = fseg.failure_mode(fy_sh=355*MPa, fy_fl=335*MPa)
print(f"Governing failure mode: {mode}")
print(f"Ultimate shell forces A, B, D, E [kN]: {[round(Z/kN) for Z in Zu]}")

fseg.validate(fy_sh=355*MPa, fy_fl=335*MPa)    # raises ValueError if mode B does not govern
# --8<-- [end:validate]


# --8<-- [start:points]
print(f"Preload after losses: {fseg.bolt_preload/kN:.0f} kN")
print(f"P1 (at rest):       Z = {fseg.shell_force_at_rest/kN:8.1f} kN   Fs = {fseg.bolt_force_at_rest/kN:7.1f} kN")
print(f"P2 (tensile ULS):   Z = {fseg.shell_force_at_tensile_ULS/kN:8.1f} kN   Fs = {fseg.bolt_force_at_tensile_ULS/kN:7.1f} kN")
print(f"P3 (small displ.):  Z = {fseg.shell_force_at_small_displacement/kN:8.1f} kN   Fs = {fseg.bolt_force_at_small_displacement/kN:7.1f} kN")
print(f"P4 (gap closed):    Z = {fseg.shell_force_at_closed_gap/kN:8.1f} kN")
# --8<-- [end:points]


# --8<-- [start:transfer]
import numpy as np

Z = np.linspace(-500*kN, 1500*kN, 201)     # shell force per bolt
Fs = fseg.bolt_axial_force(Z)              # bolt axial force
Ms = fseg.bolt_bending_moment(Z)           # bolt bending moment

print(f"Fs(Z = 500 kN) = {fseg.bolt_axial_force(500*kN)/kN:.1f} kN")
print(f"Ms(Z = 500 kN) = {fseg.bolt_bending_moment(500*kN)/1000:.2f} kNm")
# --8<-- [end:transfer]


# --8<-- [start:plot]
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None      # matplotlib is optional: py -m pip install matplotlib

if plt:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.plot(Z/kN, Fs/kN)
    ax1.set_xlabel("Shell force Z [kN]")
    ax1.set_ylabel("Bolt axial force Fs [kN]")
    ax1.grid(True)
    ax2.plot(Z/kN, Ms/1000)
    ax2.set_xlabel("Shell force Z [kN]")
    ax2.set_ylabel("Bolt bending moment Ms [kNm]")
    ax2.grid(True)
    fig.tight_layout()
    fig.savefig("lflange-transfer-functions.png")
# --8<-- [end:plot]


# --8<-- [start:fatigue]
from pyflange.fatigue import MarkovMatrix, BoltFatigueAnalysis

# Tower bending moment Markov matrix at the flange, for 30 years
mkvm = MarkovMatrix(
    range  = np.array([ 25,  50, 100, 150, 200, 250,  300, 375, 500, 625]) * 1e6,   # N·m
    mean   = np.array([100, 100, 100, 100, 100, 100,  100, 100, 100, 100]) * 1e6,   # N·m
    cycles = np.array([3e8, 1e8, 3e7, 8e6, 2e6, 5e5, 1.5e5, 4e4, 8e3, 1e3]),
    duration = 30)                                                                 # years

fatigue = BoltFatigueAnalysis(fseg, mkvm)

print(f"S-N curve knee point: {fatigue.fatigue_curve.curves[0].DS_ref/MPa:.1f} MPa")
print(f"Damage:       {fatigue.damage:.3f}")
print(f"Fatigue life: {fatigue.fatigue_life:.1f} years")
# --8<-- [end:fatigue]
