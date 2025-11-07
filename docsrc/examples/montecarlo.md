
Example: Montecarlo Simulation
============================== 

This example shows how the pythagoras library can be used to generate a random 
sample of bolt fatigue cases, which could be used, for example in a reliability 
analysis.

The first step would be to create certain 'sampler' objects, which are just
[Python generator functions](https://realpython.com/ref/glossary/generator/)
that yield random values. The samplers that we need to create are namely:

- A flange segment sampler: a random generator of `pyflange.flangesegments.FlangeSegment` objects
- A Markov matrix sampler: a random generator of flange load `pyflange.fatigue.MarkovMatrix` objects
- A SN-curve sampler: a random generator of `pyflange.fatigue.BoltFatigueCurve` objects
- An allowable damage sampler: a random generator of real numbers that represent the allowable damage

Let's start with creating a **polynomial L-Flange segment sampler**. This can be done using the `pyflange.stats.standard_PolynomialLFlangeSegment_sampler` function, as shown in the following snippet.


```py
import pyflange.stats as stats

from math import pi
mm = 0.001
kN = 1000

from pyflange.bolts import StandardMetricBolt, RoundNut
N_BOLTS = 156

fseg_samp = stats.standard_PolynomialLFlangeSegment_sampler (
        a = 150*mm,        # distance between inner face of the flange and center of the bolt hole
        b = 122*mm,        # distance between center of the bolt hole and center-line of the shell
        s =  54*mm,        # shell thickness
        t = 172*mm,        # flange thickness
        R = 4000*mm,       # shell outer curvature radius
        central_angle = 2*pi / N_BOLTS,  # angle subtended by the flange segment

        Zg = -18044*kN / N_BOLTS,      # load applied to the flange segment shell at rest

        # Bolt object representing the flange segment bolt
        bolt = StandardMetricBolt("M80", "10.9", shank_length=160*mm, 
                shank_diameter_ratio=76.1/80, stud=True),

        bolt_preload_ratio = 0.750,     # Ratio between mean preload Fp and Fy = As*fy.
        bolt_preload_cov = 0.03,        # Coefficient of variation of the bolt preload.

        Do = 86*mm,               # Bolt hole diameter
        washer = None,            # Bolt washer
        nut = RoundNut("M80"),    # Bolt nut

        flange_flatness_tolerance = 0.0014,   # 1.4 mm/m

        s_ratio = 1.0       # Ratio of bottom shell thickness over s. Default s_botom = s.
    )
```

> The `fseg_samp` object (like any other sampler) could be used in a for-loop such 
> as `for random_fseg in fseg_samp: ...` to generate an infinite (you need to break the
> loop yourself) series of random flange-segments or with the `next` function 
> (e.g. `random_fseg1 = next(fseg_samp)`, `random_fseg2 = next(fseg_samp)`, etc.) to
> generate individual random flange-segment objects. 

The next sampler we need is a **MarkovMatrix sampler**, which we can create using the
`pyflange.stats.standard_markov_matrix_sampler` fuctions, as shown in the snipped below.

```py
import pyflange.stats as stats
from pyflange.fatigue import MarkovMatrix

# An "average" flange-load markov matrix is necessary as starting point. This could be 
# created cexplicitly as shown before, but more likely it will be loaded from an external
# csv or excel file. See the pyflange.fatigue module for more details.
avg_markov_matrix = MarkovMatrix(
        cycles = [ n0,  n1,  n2, ...],
        mean   = [ M0,  M1,  M2, ...],
        range  = [DM0, DM1, DM2, ...],
        duration = 25 # years
    )

# Create a sampler that generates random markov matrices having the same number of
# cycles and mean values that the average markov matrix, but random moment ranges from a
# log-normal distribution with mean value the avg_markov_matrix range and CoV=0.12.
markov_matrix_samp = stats.standard_markov_matrix_sampler(flange_mkv)
```

> The "standard" markov matrix sampler (like any other standard sampler contained in
> the `pyflange.stats` module) is just one of the possible MarkovMatrix samplers. If you
> want you can create your own sampler and, as long as it yields MarkovMatrix objects,
> it will still work with the following code.

The third sampler we need is a **BoltFatigueCurve sampler**, which we can generate 
using the `pyflange.stats.standard_bolt_fatigue_curve_sampler` function, as shown in
the snippet below.

```py
import pyflange.stats as stats

Dn = 0.080 # Bolt nominal diameter: 80 mm
fatigue_curve_sampler = stats.bolt_fatigue_curve_sampler(Dn)
```

The forth sampler we need to create is the **allowable damage sampler**, which is 
necessay to take into account the uncertainties connected to the Miner's rule. We
can assume a sampler with log-normal distribution with mean value 1 and CoV=0.3.

```py
import pyflange.stats as stats

allowable_damage_samp = stats.lognorm_samp(1.0, 0.30)
```

The four samplers we created are actually needed to create a "master" sampler that
uses them to create random realizations of fatigue cases. This is the sampler that
we will use in the montecarlo simulation. It is created as shown in the following
snippet.

```py
import pyflange.stats as stats

fcase_samp = stats.fatigue_case_sampler(fseg_samp, 
                                        markov_matrix_samp,
                                        fatigue_curve_samp,
                                        allowable_damage_samp)
```

Each fatogue case ralization (i.e. each `fcase = next(fcase_samp)` value) will return a 
`pyflange.fatigue.BoltFatigueAnalysis` object that contains the following
attributes (see BoltFatigueAnalysis API documentation for more details):

- `fcase.fseg` : a random realization of flange segment
- `fcase.flange_mkvm` : a random realization of a flange Markov matrix
- `fcase.fatigue_curve` : a random BoltFatigueCurve realization
- `fcase.allowable_damage` : a random realization of allowable damage
- `fcase.damage` : the cumulated damage corresponding to the random attributes above
- `fcase.fatigue_life` : the fatigue life corresponding to the random attributes above

Finally, we can use the fatigue case sampler to create a sample of fatigue cases. In the
following example, we will store a 100 000 items sample in a Pandas DataFrame.
```py
import pandas as pd
SAMPLE_SIZE = int(100000)

# Create an empty DataFrame
samp_df = pd.DataFrame(dtype=float,                 # cells data type
                       index=range(SAMPLE_SIZE),    # 100000 rows
                       columns=[                    # Names of the data-frame columns
                            "BoltPreload",
                            "GapShapeFactor",
                            "TiltAngle",
                            "GapAngle",
                            "GapHeight",
                            "Damage",
                            "AllowableDamage",
                            "FatigueLife"
                       ])

# Fill the data-frame with fatigue cases, one for each row.
for i in range(SAMPLE_SIZE):

    # Generate the next random fatigue case
    fcase = next(fcase_samp)

    # Fill row number i of the data-frame with the values in fcase
    samp_df.at[i, 'BoltPreload'    ] = fcase.fseg.Fv               # Bolt preload
    samp_df.at[i, 'GapShapeFactor' ] = fcase.fseg.gap_shape_factor # Gap shape factor
    samp_df.at[i, 'TiltAngle'      ] = fcase.fseg.tilt_angle       # Flange tilt angle
    samp_df.at[i, 'GapAngle'       ] = fcase.fseg.gap_angle        # Gap angle
    samp_df.at[i, 'GapHeight'      ] = fcase.fseg.gap_height       # Gap height
    samp_df.at[i, 'Damage'         ] = fcase.damage                # Cumulate damage
    samp_df.at[i, 'AllowableDamage'] = fcase.allowable_damage      # Allowable damage
    samp_df.at[i, 'FatigueLife'    ] = fcase.fatigue_life          # Fatigue life

    # This can take a while, so print a message every 1000 iteration
    # to show the progress
    if i % 1000 == 0:
        print(f"Completed iteration number {i}")
```


The generate sample dataframe (the `samp_df` object) can be further processed with your
favourite python statictical library or exported to a .csv file or .xlsx file for
post-processing in Excel.

