# Monte Carlo from the command line

A deterministic analysis uses design values for every input. A Monte Carlo
simulation instead draws many random *realizations* of the uncertain
inputs (preload, gap, inclination, loads, S-N curve, Miner's rule), computes
the damage of each one and gives the statistical distribution of the
results, for example to estimate a probability of failure. This is the
approach used to calibrate the safety factors of IEC 61400-6 (refs. [1] and
[5], section 6).

This tutorial runs a Monte Carlo simulation with the `pyflange` command line,
without writing code.

## 1. The case file

Use the example case
[fatigue-case.json](../examples/fatigue-case.json){ download="fatigue-case.json" }
from [Getting started](../getting-started/first-analysis.md). Besides the
deterministic data, it contains a `Stats` section with the random variables:

```json
"Stats": {
    "Preload":            { "Distribution": "Normal", "MeanRatio": 0.77, "CoV": 0.03 },
    "PlasticLoss":        { "Distribution": "Normal", "Mean": 0.03, "CoV": 0.33333 },
    "DFT":                { "Distribution": "Normal", "Mean": 250e-6, "CoV": 0.12 },
    "Gap": {
        "JointDistribution": "IEC 61400-6:2020 AMD1",
        "FlatnessTolerance": 0.0014,
        "Angle":       { "Distribution": "LogNormal", "Mean": 1.745, "CoV": 1.0 },
        "ShapeFactor": { "Distribution": "Normal", "Mean": 1.0, "CoV": 0.15 }
    },
    "Tilt":               { "Distribution": "LogNormal", "Mean": 0.001745, "CoV": 0.5 },
    "LoadRangeFactor":    { "Distribution": "LogNormal", "Mean": 1.0, "CoV": 0.12 },
    "FatigueClassFactor": { "Distribution": "Normal", "Mean": 1.0, "CoV": 0.10, "DS_ref_mean": 62e6 },
    "AllowableDamage":    { "Distribution": "LogNormal", "Mean": 1.0, "CoV": 0.30 },
    "StressMultiplicationFactor": { ... }
}
```

The values follow the stochastic model of ref. [5], Table 8. Each key is
described in the [JSON format](../json-format.md#stats) reference, which
also explains [how a realization is generated](../json-format.md#monte-carlo-model).

## 2. Generate the sample

```doscon
C:\Users\yourname\flange-analyses>py -m pyflange sample fatigue-case.json --size 2000 --seed 1
Performing probabilistic analysis for case fatigue-case.json ...
... loading data from 'C:\Users\yourname\flange-analyses\fatigue-case.json'
... generating random sample of 2000 realizations with seed 1
... ... generated 1000 realizations
... ... generated 2000 realizations
... writing results to file '.\fatigue-case-sample.csv'
... DONE!
```

- `--size` is the number of realizations. 2000 takes about 15 seconds. For
  an estimate of a probability of failure of 10<sup>-4</sup> you need well
  over 10<sup>5</sup> realizations (the default is 100 000).
- `--seed` makes the simulation reproducible: the same seed and case always
  give the same sample. Without it, every run is different.
- `--xlsx` also writes the sample to an Excel file.

The CSV file has one row per realization, numbered from 0, with the random
inputs and the results (`Damage`, `AllowableDamage`, `FatigueLife`,
`Utilization` = damage / allowable damage). The columns are listed on the
[command line](../cli.md#sample) page.

## 3. Statistics

```doscon
C:\Users\yourname\flange-analyses>py -m pyflange stats fatigue-case-sample.csv
```

```
                                         MIN           P05           AVG           P95           MAX           STD       CoV
... ... BoltPreload             2.741725e+06  2.812730e+06  2.861887e+06  2.908429e+06  3.044514e+06  2.906781e+04  0.010157
... ... GapAngle                3.637145e-02  3.154774e-01  1.235215e+00  2.712764e+00  3.131882e+00  7.338239e-01  0.594086
... ... GapHeight               1.509402e-08  3.809575e-05  1.148028e-03  3.626949e-03  9.676416e-03  1.224387e-03  1.066514
... ... Damage                  1.453893e-03  8.273311e-03  1.175896e-01  3.978284e-01  1.148680e+01  3.077050e-01  2.616772
... ... Utilization             8.771545e-04  8.206894e-03  1.296999e-01  4.197863e-01  1.498221e+01  3.874562e-01  2.987328

... ... sample size: 2000
```

(Some rows are omitted here.) The values are in SI units: angles in
radians, gap heights in m, preloads in N.

To compute other statistics, open the CSV file in Excel or in Python. For
example, the fraction of realizations with `Utilization > 1` estimates the
probability of failure; in this sample it is 16 out of 2000 (0.8%).

## 4. Inspect a realization

The largest utilization of the sample, about 15, belongs to realization 928.
To see why, rebuild it and write its complete report, as `determ` does for
the deterministic case:

```doscon
C:\Users\yourname\flange-analyses>py -m pyflange inspect fatigue-case.json 928
...
... writing data to json file '.\fatigue-case-realization-928.json'
... DONE!
```

The report shows that this realization has a 92° gap of 4.8 mm, about 1.7
times the AMD1 design gap height for that angle: one of the rare large gaps
in the tail of the gap height distribution.

## Next steps

- To use different distributions, or random variables that the JSON format
  does not support, build the simulation in Python:
  [Monte Carlo in Python](montecarlo-python.md).
- To use the ED2 model and gap statistics, set
  `"type": "PolynomialLFlangeSegmentED2"` in `FlangeSegment` and
  `"JointDistribution": "IEC 61400-6:2020 ED2"` in `Stats.Gap`.
