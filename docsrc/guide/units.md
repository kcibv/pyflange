# Units and conventions

## Units

The Python API does **not** convert units: it computes with the numbers you
give it. Use one consistent system for all inputs, and the results come out
in the same system. The recommended system, used throughout this
documentation and by the standard databases, is SI:

| Quantity | Unit |
|----------|------|
| Length | m |
| Force | N |
| Moment | N·m |
| Stress, elastic modulus | Pa |
| Angle | rad |
| Stiffness of the shell (`k_shell`) | N/m per m of segment width (N/m²) |
| Duration of a Markov matrix | free (e.g. years); the fatigue life has the same unit |

!!! danger "Standard bolts and the FEA database are in SI units"
    [`StandardMetricBolt`][pyflange.bolts.StandardMetricBolt], the nuts,
    washers and hole diameters, and
    [`shell_stiffness`][pyflange.flangesegments.shell_stiffness] return
    values in m and Pa. If you use them, all other inputs must be in SI
    units too.

A convenient way to write readable input is to define unit factors:

```python
from math import pi
mm, kN, MPa, deg = 0.001, 1000, 1e6, pi/180

t = 172*mm          # 0.172 m
Fv = 2876*kN        # 2 876 000 N
gap_angle = 30*deg  # 0.5236 rad
```

The [JSON case format](../json-format.md) always uses SI units. The command
line options `--gap-angle`, `--gap-height` and `--tilt-angle` are the only
exceptions: they are given in degrees and millimeters for convenience.

## Sign conventions

- **Shell force Z**: positive in tension. The dead weight `Zg` is negative.
- **Gap closing forces** (`DZ_gap`, `DZ_gap_tot` in the reports): negative
  values close the gap.
- **Flange inclination** `tilt_angle`: positive values represent the inside
  inclination of the flange contact surface considered in ref. [2], G.2.4,
  which improves the bolt force curve (ref. [5], 6.4.2).
- **Bolt bending moment**: positive under tensile shell forces. Only its
  range matters for the fatigue analysis.

## Shell force and bending moment

The flange segment models work with the shell force **per segment**. To get
it from a tower bending moment M, PyFlange uses (ref. [2], eq. G.33):

$$Z = \frac{M}{W} \, s \, c$$

where W is the section modulus of the tower shell and s·c the shell
cross-section area of the segment. This conversion is done by
[`bolt_markov_matrix`][pyflange.flangesegments.bolt_markov_matrix] and by
[`BoltFatigueAnalysis`][pyflange.fatigue.BoltFatigueAnalysis], which expect a
Markov matrix of **tower bending moments**. To evaluate the transfer
functions directly, pass the shell force per segment.
