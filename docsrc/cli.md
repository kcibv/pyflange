Command-line interface
======================

PyFlange includes a command-line tool that performs deterministic and Monte
Carlo bolt fatigue analyses of flange cases described in the
[PyFlange JSON format](json-format.md). After installing PyFlange, it can be
invoked either as `pyflange` or as `python -m pyflange`:

```bash
pyflange --help
pyflange COMMAND --help
```

A complete example case file is available here:
[fatigue-case.json](examples/fatigue-case.json).


Commands
--------

### determ

Performs a deterministic fatigue analysis and writes a detailed report of
the case and its results to a JSON file (see [Result files](json-format.md#result-files)).

```bash
pyflange determ case.json                    # writes case-res.json
pyflange determ case.json --xlsx             # also writes case-res.xlsx
pyflange determ case.json --gap-angle 45 --gap-height 1.5 --tilt-angle 0.2 --smf 1.1
```

| Option | Description |
|--------|-------------|
| `--gap-angle` | Override of the gap angle, in **degrees**. The shell stiffness is re-interpolated for the new angle. |
| `--gap-height` | Override of the gap height, in **mm**. |
| `--tilt-angle` | Override of the tilt angle, in **degrees**. |
| `--smf` | Stress multiplication factor. Defaults to the case `SMF`, or 1.0. |
| `--odir` | Output folder. Defaults to the current folder. |
| `--oname` | Output file name. Defaults to `<case>-res[<overrides>].json`. |
| `--xlsx` | Also write the report to an Excel file, as a key/value table with a named cell per value. |

### sample

Runs a Monte Carlo simulation, based on the `Stats` section of the case, and
writes one row per realization to a CSV file.

```bash
pyflange sample case.json                     # 100000 realizations, case-sample.csv
pyflange sample case.json --size 5000 --seed 1 --xlsx
```

| Option | Description |
|--------|-------------|
| `--size` | Number of realizations. Defaults to 100000. |
| `--seed` | Seed of the random number generator. With the same seed and case, the same sample is generated. If omitted, every run is different. |
| `--odir` | Output folder. Defaults to the current folder. |
| `--oname` | Output file name. Defaults to `<case>-sample.csv`. |
| `--xlsx` | Also write the sample to an Excel file. |

The CSV columns are:

| Column | Unit | Description |
|--------|------|-------------|
| `BoltPreload` | N | Bolt preload after losses. |
| `GapShapeFactor` | - | Gap shape factor. |
| `TiltAngle` | rad | Flange tilt angle. |
| `GapAngle` | rad | Gap angle. |
| `GapHeight` | m | Gap height. |
| `FatigueClass` | Pa | Knee-point stress range of the random bolt S-N curve (size effect included). |
| `MomentRangeCoefficient` | - | Factor applied to the Markov matrix ranges. |
| `Damage` | - | Fatigue damage. |
| `AllowableDamage` | - | Allowable damage. |
| `FatigueLife` | as Markov matrix duration | Fatigue life. |
| `SMF` | - | Stress multiplication factor. |
| `Utilization` | - | `Damage / AllowableDamage`. |

### inspect

Rebuilds a single realization of a sample and writes its detailed report,
like `determ` does for the deterministic case. This is useful to investigate
outliers.

```bash
pyflange inspect case.json 1234               # reads case-sample.csv, writes case-realization-1234.json
pyflange inspect case.json 1234 --spath other-sample.csv --xlsx
```

| Option | Description |
|--------|-------------|
| `--spath` | Path of the sample CSV file. Defaults to `<case>-sample.csv`, in the folder of the case. |
| `--odir`, `--oname`, `--xlsx` | As for `determ`. |

### stats

Prints minimum, 5% quantile, average, 95% quantile, maximum, standard
deviation and coefficient of variation of each column of a sample.

```bash
pyflange stats case-sample.csv
```


Migrating from Forrester
------------------------

The PyFlange command-line interface replaces the Forrester tool of the Bolt &
Beautiful project:

- `forrester determ`, `sample`, `inspect` and `stats` become `pyflange determ`,
  `sample`, `inspect` and `stats`, with the same options. `sample` adds the
  `--seed` option and writes the same CSV columns.
- Forrester case files (format version 0.1.0) can be used as they are.
- `forrester hist` has been removed: plot the sample CSV with your tool of choice.
- `forrester export-sgre-case` is not part of PyFlange: the conversion of
  company-specific file formats to the PyFlange JSON format is left to
  dedicated tools.
