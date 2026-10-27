# First analysis without code

In this step you run a deterministic bolt fatigue analysis of an example
L-flange, using the PyFlange command-line tool. You need PyFlange
[installed](install.md).

## 1. Create a working folder

Open a terminal (++win++, type `cmd`, ++enter++) and create a folder for
your analyses:

```doscon
C:\Users\yourname>mkdir flange-analyses
C:\Users\yourname>cd flange-analyses
C:\Users\yourname\flange-analyses>
```

## 2. Download the example case

Download the example case file
[fatigue-case.json](../examples/fatigue-case.json){ download="fatigue-case.json" }
and save it in the `flange-analyses` folder. (In most browsers: right-click
the link, **Save link as...**.)

Check that the file is there:

```doscon
C:\Users\yourname\flange-analyses>dir
...
23/10/2026  10:15             3.512 fatigue-case.json
```

!!! tip "Look inside the case file"
    A case file is a plain text file. Open it with Notepad
    (`notepad fatigue-case.json`) or, better, with an editor such as
    [VS Code](https://code.visualstudio.com/) or Notepad++, which highlight
    the structure. It describes an L-flange of an 8 m tower with 156 M80
    bolts, a 30° design gap, a 30-year load spectrum and the bolt S-N curve.
    All values are in SI units (m, N, Pa, rad). Every field is explained in
    the [JSON format](../json-format.md) reference.

## 3. Run the analysis

```doscon
C:\Users\yourname\flange-analyses>py -m pyflange determ fatigue-case.json --xlsx
Performing deterministic fatigue analysis for case 'fatigue-case.json' ...
... loading data from 'fatigue-case.json'
... building the model
... writing data to json file '.\fatigue-case-res.json'
... writing data to excel file '.\fatigue-case-res.xlsx'
... DONE!
```

`determ` stands for *deterministic*. The `--xlsx` option also writes the
results to an Excel file.

## 4. Look at the results

Open the Excel file:

```doscon
C:\Users\yourname\flange-analyses>start fatigue-case-res.xlsx
```

The sheet is a two-column table of all the inputs and results. The most
important rows are:

| Key | Value | Meaning |
|-----|-------|---------|
| `damage` | 0.501 | Miner's sum of the bolt for the load spectrum. |
| `FatigueLife` | 59.8 | Fatigue life, in the same unit as the Markov matrix duration (years). |
| `FlangeSegment.BoltPreload` | 2733846 | Bolt preload after settlement and plastic losses, in N. |
| `FlangeSegment.Point1.Z` ... `Point4.Fs` | | The points of the bolt force polynomial (see [Model overview](../guide/model.md)). |
| `FlangeSegment.TransferFunctions.Z`, `.Fs`, `.Ms` | | The bolt force and moment curves, as lists of 100 values. |

Each value cell has an Excel name equal to its key, so you can refer to it
from another workbook, e.g. `='fatigue-case-res.xlsx'!damage`.

The damage is lower than the allowable damage of 1.0: the bolts of this
flange pass the fatigue check for the 30° gap.

## 5. Try another gap

IEC 61400-6 requires the gap angles 30°, 60°, 90° and 120° to be checked.
Run the same case with a 90° gap of 2.75 mm, the AMD1 design gap height for
this tower diameter:

```doscon
C:\Users\yourname\flange-analyses>py -m pyflange determ fatigue-case.json --gap-angle 90 --gap-height 2.75
```

The result is written to `fatigue-case-res-GapA90.0deg-gaph2.75mm.json`.
Open it with Notepad and look for `"damage"`: it is larger than 1, so the
90° gap governs and the bolts fail the check. The
[design gaps tutorial](../tutorials/design-gaps.md) shows how to run this
check systematically.

## What you have learned

- A flange case is described in a JSON file.
- `py -m pyflange determ CASE.json` runs the deterministic analysis, and
  the gap and tilt can be changed from the command line.
- The [command-line reference](../cli.md) lists all the commands, including
  the Monte Carlo simulation (`sample`).

## Next step

[Write your first Python script :octicons-arrow-right-24:](first-script.md)
