# PyFlange tutorial: writing and reading JSON case files from Python.
#
# This script is explained step by step in the tutorial
# https://kcibv.github.io/pyflange/latest/tutorials/case-files/
#
# Run it in a folder that contains the example case file fatigue-case.json.

from math import pi
deg = pi/180


# --8<-- [start:load]
from pyflange.io import load_inputfile, create_fatiguecase

case = load_inputfile("fatigue-case.json")      # a plain Python dictionary
fcase = create_fatiguecase(case)                # a BoltFatigueAnalysis object
print(f"Damage of the example case: {fcase.damage:.3f}")
# --8<-- [end:load]


# --8<-- [start:variants]
import copy
from pyflange.io import dict_to_json
from pyflange.stats import gap_height_distribution

D = 2 * case["FlangeSegment"]["R"]          # tower outer diameter
u_tol = 0.0014                              # flatness tolerance, 1.4 mm/m

for angle_deg in (30, 60, 90, 120):
    angle = angle_deg * deg
    variant = copy.deepcopy(case)
    variant["CaseID"] = f"Example, {angle_deg} deg gap"
    variant["FlangeSegment"]["gap"]["angle"] = angle
    variant["FlangeSegment"]["gap"]["height"] = gap_height_distribution(D, u_tol, angle * D/2).ppf(0.95)
    dict_to_json(variant, f"case-gap{angle_deg}.json")
    print(f"Written case-gap{angle_deg}.json")
# --8<-- [end:variants]


# --8<-- [start:export]
from pyflange.io import fatiguecase_to_dict

result = fatiguecase_to_dict(fcase, verbose=True)   # input + results + report
dict_to_json(result, "case-res.json")
print(f"Fatigue life: {result['FatigueLife']:.1f} years")
# --8<-- [end:export]


# --8<-- [start:batch]
import glob
import pandas as pd
from pyflange.io import load_fatiguecase

rows = []
for path in sorted(glob.glob("case-gap*.json")):
    fc = load_fatiguecase(path)
    rows.append({"case": path,
                 "gap angle [deg]": fc.fseg.gap.angle / deg,
                 "damage": fc.damage,
                 "fatigue life [years]": fc.fatigue_life})

table = pd.DataFrame(rows)
print(table.to_string(index=False))
table.to_excel("results.xlsx", index=False)
# --8<-- [end:batch]
