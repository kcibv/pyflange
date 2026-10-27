# References

The documentation and the docstrings of the source code cite the following
documents as `ref. [n]`. The numbers are the same everywhere.

| Ref. | Document | Used for |
|------|----------|----------|
| [1] | Seidel, M.: *IEC 61400-6 AMD1 Background document. Fatigue design for ring flange connections in wind turbine support structures.* Rev00, Siemens Gamesa, 2024. | Background of the AMD1 model, design gaps, stochastic model (Table 7). |
| [2] | IEC 61400-6:2020+AMD1:2025 CSV. *Wind energy generation systems - Part 6: Tower and foundation design requirements.* Consolidated version, IEC, 2025. | The standard: design preload and losses (6.7.5.1), design gap (6.7.5.2), bolt force model (6.7.5.3 and Annex G), damage calculation (6.7.5.4). Equation numbers G.xx refer to Annex G. |
| [3] | Petersen, C.: *Nachweis der Betriebsfestigkeit exzentrisch beanspruchter Ringflanschverbindungen.* Stahlbau 67 (1998), pp. 191-203. [doi:10.1002/stab.199800690](https://doi.org/10.1002/stab.199800690) | Stiffness of the clamped parts. |
| [4] | Petersen, C.: *Stahlbau.* 4th edition, Springer Vieweg, 2012. [doi:10.1007/978-3-8348-8610-1](https://doi.org/10.1007/978-3-8348-8610-1) | Stiffness of the clamped parts. |
| [5] | Seidel, M.: *IEC 61400-6 AMD1 & ED2 Background document. Fatigue design for ring flange connections in wind turbine support structures.* Rev01, draft 2026-08-04, Siemens Gamesa. | ED2 changes of the model (9.7), ED2 design gaps (2.2.4), stochastic model (6.4, Table 8), shell stiffness database (9.3.1). |
| [6] | Tobinaga, I.; Ishihara, T.: *A study of action point correction factor for L-type flanges of wind turbine towers.* Wind Energy 21 (2018), pp. 801-806. [doi:10.1002/we.2193](https://doi.org/10.1002/we.2193) | Action point correction a' (ref. [2], G.1.3). |

When a docstring cites a section or an equation, the number refers to the
document cited: for example *ref. [2], eq. G.36* is equation G.36 of
IEC 61400-6:2020+AMD1:2025, and *ref. [5], section 9.7* is section 9.7 of the
Rev01 background document.

!!! note "Availability"
    IEC 61400-6 can be purchased from the [IEC webstore](https://webstore.iec.ch)
    or from the national standardization bodies. The background documents
    are working documents of the IEC 61400-6 maintenance team; the Bolt &
    Beautiful project documentation is available from the
    [project website](https://grow-offshorewind.nl/project/bolt-beautiful).

## Product and design standards

The bolt, nut and washer models of [`pyflange.bolts`](apis/bolts.md) refer to:

- ISO 68-1:1998, *ISO general purpose screw threads - Basic profile.*
- ISO 724, *ISO general-purpose metric screw threads - Basic dimensions.*
- ISO 898-1:2013, *Mechanical properties of fasteners - Part 1: Bolts, screws and studs.*
- ISO 4032, *Hexagon regular nuts (style 1).*
- ISO 7089, *Plain washers - Normal series.*
- EN 1993-1-8:2005, *Eurocode 3 - Design of steel structures - Part 1-8: Design of joints.*
- VDI 2230 Part 1, *Systematic calculation of highly stressed bolted joints.*

## Citing PyFlange

If you use PyFlange in a report or a publication, please cite it as:

> PyFlange: Python library for the analysis of large bolted ring-flange
> connections. KCI The Engineers B.V., Siemens Gamesa Renewable Energy B.V.,
> TNO, developed within the Bolt & Beautiful GROW project.
> <https://github.com/kcibv/pyflange>, version x.y.z.

and state the version you used (`py -m pip show pyflange`).
