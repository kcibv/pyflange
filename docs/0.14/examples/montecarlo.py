# PyFlange tutorial: Monte Carlo simulation in Python.
#
# This script is explained step by step in the tutorial
# https://kcibv.github.io/pyflange/latest/tutorials/montecarlo-python/
#
# All values are in SI units: m, N, Pa, rad.

import os
from math import pi
import numpy as np

mm  = 0.001     # m
um  = 1e-6      # m
kN  = 1000      # N
MPa = 1e6       # Pa
deg = pi/180    # rad

# Number of realizations. The tutorial uses 2000 to run in a few seconds;
# use 100000 or more for a real assessment.
SAMPLE_SIZE = int(os.environ.get("PYFLANGE_SAMPLE_SIZE", 2000))


# --8<-- [start:seed]
np.random.seed(1)       # makes the simulation reproducible
# --8<-- [end:seed]


# --8<-- [start:fseg_sampler]
import pyflange.stats as stats
from pyflange.bolts import StandardMetricBolt, RoundNut, standard_hole_diameter

D = 8000*mm
n_bolts = 156
bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm,
                          shank_diameter_ratio="pitch", stud=True)

fseg_sampler = stats.standard_PolynomialLFlangeSegment_sampler(
    a = 150*mm, b = 122*mm, s = 54*mm, t = 172*mm, R = D/2,
    central_angle = 2*pi / n_bolts,
    Zg = -18044*kN / n_bolts,
    bolt = bolt,
    Do = standard_hole_diameter("M80"),
    washer = None,
    nut = RoundNut("M80"),

    # Random preload before losses: mean 2876 kN, CoV 3%
    preload_sampler = stats.norm_sampler(2876*kN, 0.03),

    # Random preload losses (ref. [5], section 6.4.1)
    ppl_sampler = stats.norm_sampler(0.03, 1/3),        # plastic loss: mean 3%, std 1%
    dft_sampler = stats.norm_sampler(250*um, 0.12),     # coating: mean 250 um, std 30 um

    # Random gap: AMD1 gap height distribution with random angle and shape
    gap_sampler = stats.standard_gap_sampler(
        flange_diameter = D,
        flange_flatness_tolerance = 0.0014,
        gap_angle_sampler = stats.lognorm_sampler(100*deg, 1.0),
        gap_shape_factor_sampler = stats.norm_sampler(1.0, 0.15)),

    # Random inside flange inclination: LN(0.1 deg, CoV 0.5)
    tilt_sampler = stats.lognorm_sampler(0.1*deg, 0.5),

    settlement_factor = 0.5)
# --8<-- [end:fseg_sampler]


# --8<-- [start:other_samplers]
from pyflange.fatigue import MarkovMatrix

design_mkvm = MarkovMatrix(
    range  = np.array([ 25,  50, 100, 150, 200, 250,  300, 375, 500, 625]) * 1e6,
    mean   = np.array([100, 100, 100, 100, 100, 100,  100, 100, 100, 100]) * 1e6,
    cycles = np.array([3e8, 1e8, 3e7, 8e6, 2e6, 5e5, 1.5e5, 4e4, 8e3, 1e3]),
    duration = 30)

# Load ranges scaled by X_FLS ~ LN(1.0, 0.12)
mkvm_sampler = stats.standard_markov_matrix_sampler(design_mkvm, 1.0, 0.12)

# Bolt S-N curve with DS_ref = X_FAT ~ N(62 MPa, 0.10)
sn_sampler = stats.standard_bolt_fatigue_curve_sampler(bolt.nominal_diameter, 1.0, 0.10, 62*MPa)

# Miner's rule uncertainty: allowable damage ~ LN(1.0, 0.30)
allowable_damage_sampler = stats.lognorm_sampler(1.0, 0.30)
# --8<-- [end:other_samplers]


# --8<-- [start:case_sampler]
fcase_sampler = stats.fatigue_case_sampler(fseg_sampler,
                                           mkvm_sampler,
                                           sn_sampler,
                                           allowable_damage_sampler)

fcase = next(fcase_sampler)         # one random BoltFatigueAnalysis
print(f"One realization: gap angle {fcase.fseg.gap.angle/deg:.0f} deg, "
      f"gap height {fcase.fseg.gap.height/mm:.2f} mm, damage {fcase.damage:.3g}")
# --8<-- [end:case_sampler]


# --8<-- [start:sample]
import pandas as pd

rows = []
for i in range(SAMPLE_SIZE):
    fcase = next(fcase_sampler)
    rows.append({
        "GapAngle":        fcase.fseg.gap.angle,
        "GapHeight":       fcase.fseg.gap.height,
        "BoltPreload":     fcase.fseg.bolt_preload,
        "Damage":          fcase.damage,
        "AllowableDamage": fcase.allowable_damage,
        "FatigueLife":     fcase.fatigue_life,
    })
    if (i+1) % 500 == 0:
        print(f"... {i+1} realizations")

sample = pd.DataFrame(rows)
# --8<-- [end:sample]


# --8<-- [start:results]
utilization = sample.Damage / sample.AllowableDamage
pf = np.mean(utilization > 1)       # fraction of failed realizations

print(sample.describe(percentiles=[0.05, 0.5, 0.95]).transpose().to_string())
print(f"Probability of failure over {design_mkvm.duration} years: {pf:.2e}")

sample.to_csv("montecarlo-sample.csv")      # for further processing, e.g. in Excel
# --8<-- [end:results]
