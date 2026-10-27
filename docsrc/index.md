---
hide:
  - navigation
---

<div class="pf-hero" markdown>

# PyFlange

<p class="pf-tagline">Open-source Python library for the fatigue analysis of large bolted
ring-flange connections, implementing the polynomial flange segment model of
IEC 61400-6 AMD1.</p>

</div>

PyFlange computes the **bolt axial force and bolt bending moment** of a
flange segment as a function of the tower shell force, including the effect of
flange **gaps** (flatness imperfections) and **flange inclination**. On top of
this model it provides the **bolt fatigue assessment** of IEC 61400-6 and the
tools for **probabilistic (Monte Carlo) analyses**.

It is meant for engineers who design or check the L-flanges and T-flanges of
wind turbine towers, monopiles and transition pieces, and who need a fast,
transparent and verifiable implementation of the standard.

```python
from pyflange.flangesegments import PolynomialLFlangeSegment

fseg = PolynomialLFlangeSegment(a=0.150, b=0.122, s=0.054, t=0.172, R=4.0, ...)
Fs = fseg.bolt_axial_force(500e3)       # bolt force for a 500 kN shell force
```

## Where to start

<div class="grid cards" markdown>

-   :material-rocket-launch-outline:{ .lg .middle } **Getting started**

    ---

    New to Python? Install PyFlange step by step on Windows and run your first
    analysis, with or without writing code.

    [:octicons-arrow-right-24: Install PyFlange](getting-started/index.md)

-   :material-school-outline:{ .lg .middle } **Tutorials**

    ---

    Worked examples: deterministic L- and T-flange analyses, design gap
    angles, AMD1 versus ED2, Monte Carlo simulations.

    [:octicons-arrow-right-24: Tutorials](tutorials/index.md)

-   :material-console:{ .lg .middle } **Command line**

    ---

    Run deterministic and Monte Carlo fatigue analyses of a JSON case file,
    without programming.

    [:octicons-arrow-right-24: Command line](cli.md)

-   :material-code-json:{ .lg .middle } **JSON case format**

    ---

    The file format that describes a flange fatigue case: geometry, bolt,
    gap, loads, S-N curve and statistics.

    [:octicons-arrow-right-24: JSON format](json-format.md)

-   :material-book-open-variant:{ .lg .middle } **User guide**

    ---

    How the model works, symbols, units and sign conventions, validation
    against reference results, troubleshooting.

    [:octicons-arrow-right-24: Model overview](guide/model.md)

-   :material-api:{ .lg .middle } **API reference**

    ---

    Every module, class and function, generated from the source code, with
    references to the IEC equations.

    [:octicons-arrow-right-24: API reference](apis/index.md)

</div>

## What PyFlange can do

| Task | Where |
|------|-------|
| Bolt force and bolt moment curves of L-flanges and T-flanges (IEC 61400-6 AMD1, Annex G.2) | [`flangesegments`](apis/flangesegments.md), [L-flange tutorial](tutorials/lflange.md), [T-flange tutorial](tutorials/tflange.md) |
| The changes of the model proposed for IEC 61400-6 ED2 | [Design gaps and ED2 tutorial](tutorials/design-gaps.md) |
| Check that the polynomial model is applicable (governing ULS failure mode) | [`failure_mode`, `validate`](apis/flangesegments.md) |
| Design gap heights (log-normal AMD1, Weibull ED2) | [`stats.gap_height_distribution`](apis/stats.md) |
| Shell stiffness from a database of 2160 FEA results | [`flangesegments.shell_stiffness`](apis/flangesegments.md) |
| Preload losses due to settlement and plastic strain | [`PolynomialLFlangeSegment`](apis/flangesegments.md) |
| Bolt S-N curve, size effect, Miner's rule, fatigue life | [`fatigue`](apis/fatigue.md) |
| Standard metric bolts, nuts and washers | [`bolts`](apis/bolts.md) |
| Monte Carlo simulation of the bolt fatigue damage | [Monte Carlo tutorials](tutorials/montecarlo-cli.md), [`stats`](apis/stats.md) |
| Case files, batch runs, results in Excel | [Command line](cli.md), [`io`](apis/io.md) |

## Background

### The Bolt & Beautiful project

PyFlange has been developed within **Bolt & Beautiful**, a
[GROW](https://grow-offshorewind.nl/project/bolt-beautiful) research project
on the bolted connections of offshore wind turbine support structures, by
**KCI**, **Siemens Gamesa**, **JBO** and **TNO**. The project is funded with
subsidy from the Top Sector Energy of the Dutch Ministry of Economic Affairs.

Within the project, extensive finite element studies were carried out. They
supported the calibration of the safety factors and small improvements of
the analytical flange model, which are proposed for the next edition of
IEC 61400-6 (ref. [5], section 1.1).

### The link with IEC 61400-6

IEC 61400-6 is the international standard for the design of wind turbine
towers and foundations. Its amendment 1 (**AMD1**, published in June 2025)
replaced the former Schmidt/Neuper approach for the fatigue design of flange
bolts with a new analytical model, developed by Marc Seidel (Siemens Gamesa),
which covers the effect of flange gaps of the size allowed by the
manufacturing tolerances (ref. [1]). The model is given in Annex G.2 of the
standard (ref. [2]), and the use of AMD1 becomes mandatory for IEC
certification from December 2026 (ref. [5], section 9.2).

The model is considerably more complex than the previous one. To help the
industry implement it, PyFlange was developed as an open-source reference
implementation, and it is referred to in the background document of the
standard (ref. [5], section 9.2). The results of PyFlange are verified against
the reference spreadsheets of Siemens Gamesa (see
[Validation](guide/validation.md)).

The revision of the standard towards Edition 2 (**ED2**) introduces some small
changes to the model and new design gap heights (ref. [5], section 9.7).
PyFlange implements them in separate classes (`PolynomialLFlangeSegmentED2`,
`PolynomialTFlangeSegmentED2`, `gap_height_distribution_ED2`), so that both
versions can be compared.

!!! warning "ED2 is a draft"
    The ED2 implementation follows the draft background document of
    2026-08-04 (ref. [5]). It may change before ED2 is published. Use the AMD1
    classes for certified designs.

The documents cited as `ref. [n]` throughout this documentation are listed on
the [References](references.md) page.

## License and disclaimer

PyFlange is free software, released under the
[GNU General Public License v3](license.md). It is provided without any
warranty: the user remains responsible for the verification of the results
and for their use in design.

Questions, bug reports and suggestions are welcome as
[GitHub issues](https://github.com/kcibv/pyflange/issues). See
[Contributing](contributing.md).
