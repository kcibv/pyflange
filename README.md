> This software is part of the Bolt and Beautiful project funded with subsidy
> from the Top Sector Energy of the Dutch Ministry of Economic Affairs.


PyFlange
=========================================================================

PyFlange is an open-source Python library for the fatigue analysis of large
bolted ring-flange connections, such as the L-flanges and T-flanges of
offshore wind turbine towers, monopiles and transition pieces.

It implements the polynomial flange segment model of **IEC 61400-6:2020+AMD1:2025**
(Annex G.2, developed by Marc Seidel), which predicts the bolt axial force and
bolt bending moment as functions of the tower shell force, including the
effect of flange gaps and flange inclination. On top of it, PyFlange provides:

- the L-flange and T-flange models of AMD1, and the changes proposed for ED2;
- the check of the model validity (governing ULS failure mode);
- standard metric bolts, nuts and washers;
- design gap height distributions and a database of FEA shell stiffnesses;
- the bolt fatigue assessment of IEC 61400-6 (S-N curve, Miner's rule);
- random samplers for Monte Carlo simulations of the bolt fatigue damage;
- a JSON case format and a command-line tool (`pyflange`) to run
  deterministic and Monte Carlo analyses without writing code.

PyFlange was developed within the Bolt & Beautiful GROW project by KCI,
Siemens Gamesa, JBO and TNO.

**Documentation: <https://kcibv.github.io/pyflange/>**


Installation
-------------------------------------------------------------------------

PyFlange requires Python 3.9 or later:

```
py -m pip install pyflange
```

(`python -m pip install pyflange` on Linux and macOS.) New to Python? The
[Getting started](https://kcibv.github.io/pyflange/latest/getting-started/) guide
explains every step on Windows.


Quick start
-------------------------------------------------------------------------

### Command line

Download the example case
[fatigue-case.json](docsrc/examples/fatigue-case.json) and run:

```
py -m pyflange determ fatigue-case.json --xlsx         # deterministic fatigue analysis
py -m pyflange sample fatigue-case.json --seed 1       # Monte Carlo simulation
py -m pyflange stats fatigue-case-sample.csv           # statistics of the sample
```

See the [command-line documentation](https://kcibv.github.io/pyflange/latest/cli/)
and the [JSON case format](https://kcibv.github.io/pyflange/latest/json-format/).

### Python

```python
from math import pi
from pyflange.bolts import StandardMetricBolt, RoundNut
from pyflange.flangesegments import PolynomialLFlangeSegment, Gap, shell_stiffness
from pyflange.stats import gap_height_distribution

mm, kN, MPa, deg = 0.001, 1000, 1e6, pi/180     # SI units throughout

D = 8000*mm                                     # tower outer diameter
gap_angle = 30*deg
gap_height = gap_height_distribution(D, 0.0014, gap_angle*D/2).ppf(0.95)   # AMD1 design gap

fseg = PolynomialLFlangeSegment(
    a = 150*mm, b = 122*mm, s = 54*mm, t = 172*mm, R = D/2,
    central_angle = 2*pi/156,                   # 156 bolts
    Zg = -18044*kN / 156,                       # dead weight per bolt
    bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm, stud=True),
    Fv = 2876*kN,                               # preload before losses
    Do = 86*mm, washer = None, nut = RoundNut("M80"),
    gap = Gap(height=gap_height, angle=gap_angle),
    k_shell = shell_stiffness(D/2, 54*mm, gap_angle))

fseg.validate(fy_sh=355*MPa, fy_fl=335*MPa)     # failure mode B must govern

Fs = fseg.bolt_axial_force(500*kN)              # bolt force for a 500 kN shell force
Ms = fseg.bolt_bending_moment(500*kN)           # bolt bending moment
```

The fatigue analysis takes a Markov matrix of tower bending moments:

```python
import numpy as np
from pyflange.fatigue import MarkovMatrix, BoltFatigueAnalysis

mkvm = MarkovMatrix(range    = np.array([100e6, 300e6, 500e6]),   # N·m
                    mean     = np.array([100e6, 100e6, 100e6]),   # N·m
                    cycles   = np.array([3e7, 1.5e5, 8e3]),
                    duration = 30)                                # years

fatigue = BoltFatigueAnalysis(fseg, mkvm)
print(fatigue.damage, fatigue.fatigue_life)
```

The [tutorials](https://kcibv.github.io/pyflange/latest/tutorials/) cover
L-flanges, T-flanges, design gap angles, ED2 and Monte Carlo simulations.


Development
-------------------------------------------------------------------------

```
git clone https://github.com/kcibv/pyflange.git
cd pyflange
py -m pip install -e . pytest -r docsrc/requirements.txt

py -m pytest                                                    # unit tests and tutorial examples
py validations/validate.py validations/references/*.json        # validation against the SGRE references
mkdocs serve                                                    # documentation preview
```

The documentation source is in `docsrc/`. GitHub workflows run the tests at
every push and publish the documentation, one version per release, to the
`docs/` folder of `main`, served by GitHub Pages. See the
[developer guide](https://kcibv.github.io/pyflange/latest/contributing/#developer-guide).


Contributing
-------------------------------------------------------------------------
You can contribute to this project by reporting a bug, highlighting a
necessary improvement, submitting a code improvement or by asking a
question. For instructions about how to do all these things, please
see our [contribution guidelines](CONTRIBUTING.md).

> Since the typical PyFlange user is expected not to be a professional
> programmer and probably not familiar with Git, we have created a
> [Git tutorial for non-programmers](https://kcibv.github.io/git-tutorial/)
> to make contributing a less intimidating process.


License
-------------------------------------------------------------------------
pyFlange - python library for large flanges design
Copyright (C) 2024  KCI The Engineers B.V.,
                    Siemens Gamesa Renewable Energy B.V.,
                    Nederlandse Organisatie voor toegepast-natuurwetenschappelijk onderzoek TNO.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License, as published by
the Free Software Foundation, either version 3 of the License, or any
later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License version 3 for more details.

You should have received a copy of the GNU General Public License
version 3 along with this program.  If not, see <https://www.gnu.org/licenses/>.
