Command-line interface
======================

PyFlange includes a command-line tool that performs deterministic and Monte
Carlo bolt fatigue analyses of flange cases described in the
[PyFlange JSON format](json-format.md), without writing any code. After
[installing PyFlange](getting-started/install.md), it can be invoked either
as `pyflange` or as `py -m pyflange` (`python -m pyflange` on Linux and
macOS):

```doscon
py -m pyflange --help
py -m pyflange COMMAND --help
py -m pyflange --version
```

!!! tip
    `py -m pyflange` always works. The shorter `pyflange` works only if the
    folder of the Python scripts is in your `PATH`. The examples below use
    the short form.

A complete example case file is available here:
[fatigue-case.json](examples/fatigue-case.json){ download="fatigue-case.json" }.
For a guided introduction, see
[First analysis without code](getting-started/first-analysis.md) and
[Monte Carlo from the command line](tutorials/montecarlo-cli.md).

| Command | Purpose |
|---------|---------|
| [`determ`](#determ) | Deterministic fatigue analysis of a case. |
| [`sample`](#sample) | Monte Carlo simulation of a case. |
| [`inspect`](#inspect) | Detailed report of one realization of a Monte Carlo sample. |
| [`stats`](#stats) | Summary statistics of a Monte Carlo sample. |


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

The main results in the report are `damage`, `FatigueLife` and, in the
`FlangeSegment` section, `BoltPreload`, the polynomial points `Point1` to
`Point4` and the curves `TransferFunctions`.

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

The first column (without header) is the number of the realization, starting
from 0, to be used with `inspect`.

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


Running batches
---------------

To analyze all the case files of a folder, use a loop of the Windows command
prompt:

```doscon
for %f in (*.json) do pyflange determ "%f" --odir results
```

(In a `.bat` file, write `%%f` instead of `%f`.) The `results` folder must
exist; create it with `mkdir results`. Note that the loop also picks up the
result files if they are written in the same folder as the cases.

In PowerShell:

```powershell
Get-ChildItem *.json | ForEach-Object { pyflange determ $_.FullName --odir results }
```

For more complex batches, for example to collect the results in one table,
use Python: see [Writing case files from Python](tutorials/case-files.md).


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
