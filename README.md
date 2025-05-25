> This software is part of the Bolt and Beautiful project funded with subsidy 
> from the Top Sector Energy of the Dutch Ministry of Economic Affairs.


PyFlange Python Package
=========================================================================

This package has the ambitious goal of providing all the tools engineers 
need for the design of large bolted flanges such as the flanges used in
offshore wind for connecting the turbine tower to the foundation.

Far from achieving its goal, this package currently contains only an
implementation of Marc Seidel's polynomial model for predicting bolt
forces and moments due to the tower shell force, (soon to be) published 
in the the standard IEC 61400-6-AMD1.

This package has beend developed within the Bolt and Beoutiful GROW
project by KCI, Siemens Gamesa and TNO.

The rest of this documentation will show how to get started and where to        
find extra documentation.


Getting Started
-------------------------------------------------------------------------
You can install PyFlange via pip as follows:

```
pip install pyflange
```

After installing the package, you can import it in your python code and start
using it. First of all, you need to create a `FlangeSegment` object as shown
below.

``` python

# Create the bolt object
from pyflange.bolts import StandardMetricBolt, ISOFlatWasher, ISOHexNut
M80_bolt   = StandardMetricBolt("M80", "10.9", shank_length=0.270, stud=True)
M80_washer = ISOFlatWasher("M80")
M80_nut    = ISOHexNut("M80")

# Define the gap parameters
from pyflange.gap import gap_height_distribution
from math import pi
D = 7.50                        # meters, flange outer diameter
gap_angle = pi/6                # 30 deg gap angle
gap_length = gap_angle * D/2    # outer length of the gap
u_tol = 0.0014                  # flatness tolerance in mm/mm
gap_dist = gap_height_distribution(D, u_tol, gap_length)    # lognormal distribution

# Create the FlangeSegment model
from pyflange.flangesegments import PolynomialLFlangeSegment
Nb = 120    # number of bolts
fseg = PolynomialLFlangeSegment(
    a = 0.2325,              # distance between inner face of the flange and center of the bolt hole
    b = 0.1665,              # distance between center of the bolt hole and center-line of the shell
    s = 0.0720,              # shell thickness
    t = 0.2000,              # flange thickness
    R = D/2,                 # shell outer curvature radius
    central_angle = 2*pi/Nb, # angle subtented by the flange segment arc

    Zg = -14795000 / Nb,     # load applied to the flange segment shell at rest
                             # (normally dead weight of tower + RNA, divided by the number of bolts)

    bolt = M80_bolt,         # bolt object created above
    Fv = 2876000,            # design bolt preload, after preload losses

    Do = 0.086,              # bolt hole diameter
    washer = M80_washer,     # washer object created above
    nut = M80_nut,           # nut object created above

    gap_height = gap_dist.ppf(0.95),    # maximum longitudinal gap height, 95% quantile
    gap_angle = gap_angle)              # longitudinal gap length

# Assert if the flange-segment fails with failure mode B.
# If not, an exception will be raised. 
fseg.validate(470e6, 450e6)
```

> Notice that a consistent set of units of measurements has been used for inputs, namely:
> meter for distances, radians for angles and newton for forces. It is not required to
> always use these units (meter, newton), but you should choose your units and always
> apply them consistently.

Once you have your `fseg` object, you can obtain the bolt forces and moments as follows:

``` python
Fs = fseg.bolt_axial_force(3500)    # bolt force corresponding to the tower shell force Z = 3500 N
Ms = fseg.bolt_bending_moment(2000) # bolt bending moment corresponding to the tower shell force Z = 2000 N
```

The argumment `Z`, passed to `bolt_axial_force` and `bolt_bending_moment` can also be a
numpy array. In that case an array of Fs and Ms values will be returned.

``` python
import numpy as np
Z = np.array([2000, 2500, 3000])
Fs = fseg.bolt_axial_force(Z)       # return the numpy array (Fs(2000), Fs(2500), Fs(3000))
Ms = fseg.bolt_bending_moment(Z)    # return the numpy array (Ms(2000), Ms(2500), Ms(3000))
```



Learn More
-------------------------------------------------------------------------

For more details, read the [pyflage API documentation](https://kcibv.github.io/pyflange/). 



Testing
-------------------------------------------------------------------------

Once you clone this package, you can run the unit tests (assuming you
have already the pytest module installed) as follows:

``` python
cd <path-to-package>
py -m pytest
```


Contributing
-------------------------------------------------------------------------
You can contribute to this project by reporting a bug, highlighting a
necessary improvement, submitting a code improvement or by asking a
question. For instructions about how to do all these things, please
see our [contribution guidelines](CONTRIBUTING.md).

> Since the typical PyFlange user  is expected not to be a professional
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


