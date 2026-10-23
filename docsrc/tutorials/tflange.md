# T-flange, deterministic analysis

This tutorial repeats the [L-flange tutorial](lflange.md) for a T-flange and
highlights the differences. Read the L-flange tutorial first.

**Download the full script:**
[tflange.py](../examples/tflange.py){ download="tflange.py" }

## Differences with the L-flange

| | L-flange | T-flange |
|-|----------|----------|
| Class | [`PolynomialLFlangeSegment`][pyflange.flangesegments.PolynomialLFlangeSegment] | [`PolynomialTFlangeSegment`][pyflange.flangesegments.PolynomialTFlangeSegment] |
| Bolts per segment | 1 | 2 (one on each side of the shell) |
| `central_angle` | 2π / number of bolts | 2π / number of bolt **pairs** |
| `Zg` | Dead weight per bolt | Dead weight per bolt **pair** |
| `a` | Inner face of the flange to bolt axis | Edge of the flange to bolt axis |
| `tilt_angle` | Optional (default 0) | **Required**: IEC 61400-6 requires the flange taper to be always considered (ref. [2], G.2.4) |
| Governing failure mode for the model to be valid | B | A (ref. [2], eq. G.31) |
| Tensile bolt moment | Quadratic in Z | Cubic in Z (ref. [2], eq. G.61) |
| JSON case format and command line | Supported | Not yet supported |

## 1. Flange segment

<figure>
--8<-- "assets/img/tflange-geometry.svg"
</figure>

```python
--8<-- "examples/tflange.py:segment"
```

The bolt force `Fs` and moment `Ms` returned by the model are those of
**one** bolt of the pair, for the shell force `Z` acting on the whole
segment.

`tilt_angle = 0.0` is allowed only if the taper of the flange ensures that
the first contact point is directly below the tower wall (ref. [2],
6.7.5.3). Otherwise, use the worst-case inclination allowed by the
manufacturing specification.

## 2. Validity check

```python
--8<-- "examples/tflange.py:validate"
```

```
Governing failure mode: A
Ultimate shell forces A, B, D, E [kN]: [2121, 2950, 4843, 9790]
```

For T-flanges the polynomial model requires bolt failure (mode A) to govern.

## 3. Bolt force and moment

```python
--8<-- "examples/tflange.py:transfer"
```

```
Z =      0 kN   Fs =  935.5 kN   Ms =      4 Nm
Z =    250 kN   Fs =  956.4 kN   Ms =     49 Nm
Z =    500 kN   Fs =  981.9 kN   Ms =    112 Nm
Z =   1000 kN   Fs = 1046.8 kN   Ms =    284 Nm
```

## 4. Fatigue analysis

The fatigue analysis is the same as for the L-flange:

```python
--8<-- "examples/tflange.py:fatigue"
```

```
Damage:       0.667
Fatigue life: 38 years
```

## ED2

[`PolynomialTFlangeSegmentED2`][pyflange.flangesegments.PolynomialTFlangeSegmentED2]
has the same parameters and implements the ED2 gap stiffness (ref. [5],
9.7, change #2). The other ED2 changes concern L-flanges only.
