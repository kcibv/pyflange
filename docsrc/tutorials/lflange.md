# L-flange, deterministic analysis

This tutorial performs the deterministic bolt fatigue assessment of an
L-flange according to IEC 61400-6 AMD1, for a single design gap. It covers
the whole chain: bolt, gap, flange segment, validity check, bolt force and
moment curves, fatigue damage.

**Download the full script:**
[lflange.py](../examples/lflange.py){ download="lflange.py" }. Run it with
`py lflange.py`.

## 1. Units

PyFlange works with any consistent system of units. We use SI units and
define a few factors to make the input readable:

```python
--8<-- "examples/lflange.py:units"
```

## 2. Bolt, nut and bolt hole

```python
--8<-- "examples/lflange.py:bolt"
```

- [`StandardMetricBolt`][pyflange.bolts.StandardMetricBolt] creates a bolt
  with ISO thread dimensions and the properties of the material grade.
  `stud=True` models a stud bolt with two nuts; `shank_diameter_ratio="pitch"`
  gives a reduced shank with the diameter of the thread pitch diameter.
- [`RoundNut`][pyflange.bolts.RoundNut] is a flanged round nut. Use
  [`ISOHexNut`][pyflange.bolts.ISOHexNut] for a hexagon nut, and pass a
  [`ISOFlatWasher`][pyflange.bolts.ISOFlatWasher] as `washer` if there is
  one.
- [`standard_hole_diameter`][pyflange.bolts.standard_hole_diameter] returns
  the standard clearance hole.

## 3. Design gap

The polynomial model accounts for a gap between the two flanges, caused by
flatness imperfections. IEC 61400-6 AMD1 defines the design gap height as the
95% quantile of a log-normal distribution that depends on the tower
diameter, the flatness tolerance and the gap length (ref. [2], 6.7.5.2):

```python
--8<-- "examples/lflange.py:gap"
```

For an 8 m tower with a 1.4 mm/m tolerance and a 30° gap, the design gap
height is 0.89 mm. In this tutorial we check the 30° gap only. The standard
requires the gap angles 30°, 60°, 90° and 120°: see
[Design gap angles and ED2](design-gaps.md).

## 4. Flange segment

The model describes one *segment* of the flange: the slice that contains one
bolt (see [Model overview](../guide/model.md) and
[Nomenclature](../guide/nomenclature.md) for the symbols).

<figure>
--8<-- "assets/img/lflange-geometry.svg"
</figure>

```python
--8<-- "examples/lflange.py:segment"
```

Some parameters deserve attention:

- **`Zg`** is the shell force per bolt in the unloaded state, normally the
  dead weight above the flange divided by the number of bolts. It is
  negative because it compresses the flange.
- **`Fv`** is the preload *before* losses. The losses due to settlement
  (`DFT`, `settlement_factor`, ref. [2], eqs. (27) and (28)) and plastic
  strain (`PPL`, ref. [2], eq. (24)) are subtracted by the model. The
  resulting preload is available as `fseg.bolt_preload`. If your `Fv` already
  includes all losses, leave `PPL` and `DFT` at their default value of zero.
- **`k_shell`** is the shell stiffness. Here it is interpolated from a
  database of FEA results with
  [`shell_stiffness`][pyflange.flangesegments.shell_stiffness]. If omitted,
  the simplified formula of ref. [2], eq. G.43 is used, which is
  conservative.
- **`tilt_angle`** is the inclination of the flange contact surface. A
  positive value (inclination towards the inside) reduces the bolt forces.
- `E`, `G`, `s_ratio` (bottom shell thickness over `s`) and `r` (fillet
  radius) have default values; see
  [`PolynomialLFlangeSegment`][pyflange.flangesegments.PolynomialLFlangeSegment].

## 5. Check that the model is applicable

The polynomial model of an L-flange is valid only if failure mode B (plastic
hinges in the flange and in the shell) is the governing ultimate failure
mode (ref. [2], eq. G.30):

```python
--8<-- "examples/lflange.py:validate"
```

```
Governing failure mode: B
Ultimate shell forces A, B, D, E [kN]: [3128, 1805, 2206, 5503]
```

`validate()` raises a `ValueError` if another mode governs, so a script
stops before producing results outside the validity range of the model.

## 6. Bolt force curve

The bolt force curve is defined by four points (ref. [2], Table G.1):

```python
--8<-- "examples/lflange.py:points"
```

```
Preload after losses: 2669 kN
P1 (at rest):       Z =   -115.7 kN   Fs =  2669.2 kN
P2 (tensile ULS):   Z =   1840.0 kN   Fs =  3336.6 kN
P3 (small displ.):  Z =    111.1 kN   Fs =  2677.0 kN
P4 (gap closed):    Z =   -427.5 kN
```

The transfer functions accept a single value or a numpy array of shell
forces:

```python
--8<-- "examples/lflange.py:transfer"
```

With matplotlib installed, the script also plots the curves:

```python
--8<-- "examples/lflange.py:plot"
```

![Bolt force and bolt moment curves](../assets/img/lflange-transfer-functions-light.png#only-light)
![Bolt force and bolt moment curves](../assets/img/lflange-transfer-functions-dark.png#only-dark)

Below P1 the flange is compressed and the gap closes: at P4 it is completely
closed and the bolt force stays constant. Above P1 the bolt force increases
quadratically up to the ultimate state P2. P3 sets the initial slope of the
curve.

## 7. Fatigue analysis

The fatigue loads are given as a Markov matrix of **tower bending moments**
at the flange. [`BoltFatigueAnalysis`][pyflange.fatigue.BoltFatigueAnalysis]
converts them into bolt stresses through the transfer functions and computes
the damage with the IEC bolt S-N curve (detail category 50, size effect,
γ<sub>M</sub> = 1.1, ref. [2], 6.7.5.4):

```python
--8<-- "examples/lflange.py:fatigue"
```

```
S-N curve knee point: 40.1 MPa
Damage:       0.558
Fatigue life: 53.8 years
```

The damage is below 1.0, so the bolts pass for the 30° gap. The bolt stress
Markov matrix is available as `fatigue.bolt_mkvm`, and a different S-N curve
can be passed as `custom_fatigue_curve`. A stress multiplication factor
(`SMF`) can be applied to the bolt stress ranges.

## Next steps

- Check all the design gap angles:
  [Design gap angles and ED2](design-gaps.md).
- Save this case to a JSON file and run it from the command line:
  [Writing case files from Python](case-files.md).
- The meaning of each intermediate quantity: [Model overview](../guide/model.md)
  and the [API reference](../apis/flangesegments.md).
