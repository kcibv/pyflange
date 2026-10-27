# Design gap angles and ED2

IEC 61400-6 does not prescribe a single design gap: the gap angles 30°, 60°,
90° and 120° are potentially design-driving, and the one giving the highest
fatigue damage is decisive (ref. [2], 6.7.5.3). The design gap height
depends on the gap angle. This tutorial runs the four cases for the flange of
the [L-flange tutorial](lflange.md), first with the AMD1 model and then with
the model proposed for ED2.

**Download the full script:**
[design_gaps.py](../examples/design_gaps.py){ download="design_gaps.py" }

## 1. Common parameters

All the parameters except the gap and the shell stiffness are the same for
the four cases. We store them in a dictionary and pass it to the class with
`**common`:

```python
--8<-- "examples/design_gaps.py:common"
```

## 2. AMD1

For each angle, the design gap height is the 95% quantile of the AMD1
log-normal distribution, and the shell stiffness is interpolated for that
angle:

```python
--8<-- "examples/design_gaps.py:amd1"
```

```
AMD1: design gap = 95% quantile of the log-normal distribution
  gap angle  30 deg   gap height 0.89 mm   damage 0.558
  gap angle  60 deg   gap height 1.68 mm   damage 0.869
  gap angle  90 deg   gap height 2.75 mm   damage 1.337
  gap angle 120 deg   gap height 4.12 mm   damage 1.231
  governing gap angle: 90 deg
```

The design gap heights match Table 2 of ref. [1]. The 30° gap of the
L-flange tutorial is not the governing one: with a 90° gap the damage
exceeds 1.0 and the flange fails the check.

## 3. ED2

The model proposed for ED2 (ref. [5], 9.7) changes the gap stiffness, the
gap closing force (torsional correction) and the lower bound of the
parameter a\*. The ED2 design gap height is the **75%** quantile of a
**Weibull** distribution (ref. [5], 2.2.4):

```python
--8<-- "examples/design_gaps.py:ed2"
```

```
ED2: design gap = 75% quantile of the Weibull distribution
  gap angle  30 deg   gap height 0.57 mm   damage 0.418
  gap angle  60 deg   gap height 1.27 mm   damage 0.569
  gap angle  90 deg   gap height 2.16 mm   damage 0.760
  gap angle 120 deg   gap height 3.25 mm   damage 0.736
```

With the smaller ED2 design gaps the same flange passes. Note that the ED2
design gap fractile was selected together with the partial safety factors
(ref. [5], 6.7), so a lower damage does not simply mean that ED2 is less
conservative.

!!! warning
    The ED2 classes follow the draft background document of 2026-08-04
    (ref. [5]) and may change before ED2 is published.

## Running the same study from the command line

The `determ` command accepts `--gap-angle` and `--gap-height`:

```doscon
py -m pyflange determ case.json --gap-angle 90 --gap-height 2.75
```

When the gap angle is changed, the shell stiffness is interpolated again
for the new angle. See the [command line](../cli.md) page and
[Writing case files from Python](case-files.md) for batches.
