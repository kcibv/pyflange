# Tutorials

Each tutorial is a complete, working example that you can download and run.
The scripts are tested automatically against the current version of
PyFlange, so the code on these pages always works.

| Tutorial | You will learn to | Needs |
|----------|-------------------|-------|
| [L-flange, deterministic](lflange.md) | Build an L-flange segment, check the validity of the model, plot the bolt force and moment curves, compute the fatigue damage. | Python |
| [T-flange, deterministic](tflange.md) | The same for a T-flange, and the differences with L-flanges. | Python |
| [Design gap angles and ED2](design-gaps.md) | Check the four design gap angles of IEC 61400-6 and compare the AMD1 and ED2 models. | Python |
| [Monte Carlo from the command line](montecarlo-cli.md) | Run a probabilistic analysis of a JSON case, inspect the statistics and single realizations. | Command line |
| [Monte Carlo in Python](montecarlo-python.md) | Build a custom Monte Carlo simulation with the random samplers of `pyflange.stats`. | Python |
| [Writing case files from Python](case-files.md) | Create case variants, run batches and collect the results in Excel. | Python |

If you have never used Python, start with
[Getting started](../getting-started/index.md).

## Conventions used in the tutorials

- All values are in SI units (m, N, Pa, rad). The scripts define `mm`, `kN`,
  `MPa` and `deg` to make the input readable: `150*mm` is 0.150 m.
- The example flange is an L-flange of an 8 m monopile or tower with 156 M80
  bolts. The values are realistic but do not belong to a real project.
- The Markov matrices are invented for the examples. In a real assessment,
  they come from the load calculation of the turbine.
- References such as *ref. [2], eq. G.36* point to the documents on the
  [References](../references.md) page.

!!! tip "Plots"
    Some tutorials draw plots with [matplotlib](https://matplotlib.org/),
    which is not installed together with PyFlange. Install it with
    `py -m pip install matplotlib`. The scripts work without it, they only
    skip the plots.
