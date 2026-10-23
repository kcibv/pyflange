# Monte Carlo in Python

This tutorial builds a Monte Carlo simulation of the bolt fatigue damage with
the random samplers of [`pyflange.stats`](../apis/stats.md). It does the same
as the [`sample` command](montecarlo-cli.md), but in a script, so you can
change any distribution, add random variables, or post-process the results
as you like.

**Download the full script:**
[montecarlo.py](../examples/montecarlo.py){ download="montecarlo.py" }

## Samplers

A *sampler* is a Python
[generator](https://docs.python.org/3/glossary.html#term-generator): an
object that returns a new random value every time you call `next()` on it.

```python
from pyflange.stats import norm_sampler

preload = norm_sampler(2876e3, 0.03)    # normal, mean 2876 kN, CoV 3%
next(preload)                           # a random preload, e.g. 2.91e6
next(preload)                           # another one, e.g. 2.84e6
```

PyFlange has generic samplers (`norm_sampler`, `lognorm_sampler`,
`weibull_sampler`, `beta_sampler`, `sampler` for any scipy distribution) and
samplers of PyFlange objects, which combine other samplers:

| Sampler | Yields |
|---------|--------|
| [`standard_gap_sampler`][pyflange.stats.standard_gap_sampler] | Random `Gap` (AMD1 gap height distribution) |
| [`standard_gap_sampler_ED2`][pyflange.stats.standard_gap_sampler_ED2] | Random `Gap` (ED2 distributions) |
| [`standard_PolynomialLFlangeSegment_sampler`][pyflange.stats.standard_PolynomialLFlangeSegment_sampler] | Random L-flange segments |
| [`standard_markov_matrix_sampler`][pyflange.stats.standard_markov_matrix_sampler] | Markov matrices with random load ranges |
| [`standard_bolt_fatigue_curve_sampler`][pyflange.stats.standard_bolt_fatigue_curve_sampler] | Random bolt S-N curves |
| [`fatigue_case_sampler`][pyflange.stats.fatigue_case_sampler] | Random `BoltFatigueAnalysis` objects |

You can write your own samplers too: any generator that yields objects of the
right type works.

## 1. Reproducibility

The samplers use the global numpy random generator. Seed it to get the same
sample at every run:

```python
--8<-- "examples/montecarlo.py:seed"
```

## 2. Flange segment sampler

```python
--8<-- "examples/montecarlo.py:fseg_sampler"
```

For each realization, the sampler draws a gap and a flange inclination, then
draws the preload and its losses for half of the bolts within the gap and
averages them (ref. [5], 4.4). The yielded segment has the averaged preload
after losses as `Fv`, and `PPL = DFT = 0`.

The distributions are those of the stochastic model of ref. [5], Table 8,
except the gap angle, which the reference does not give.

## 3. Loads, S-N curve and Miner's rule

```python
--8<-- "examples/montecarlo.py:other_samplers"
```

## 4. Fatigue case sampler

The four samplers are combined into a sampler of complete fatigue analyses:

```python
--8<-- "examples/montecarlo.py:case_sampler"
```

```
One realization: gap angle 46 deg, gap height 0.19 mm, damage 0.0397
```

Each realization is a [`BoltFatigueAnalysis`][pyflange.fatigue.BoltFatigueAnalysis]:
`fcase.fseg` is the random segment, `fcase.flange_mkvm` the random Markov
matrix, `fcase.fatigue_curve` the random S-N curve,
`fcase.allowable_damage` the random allowable damage, and `fcase.damage`
and `fcase.fatigue_life` the results.

## 5. Generate the sample

Collect the values of interest in a [pandas](https://pandas.pydata.org/)
table:

```python
--8<-- "examples/montecarlo.py:sample"
```

## 6. Results

```python
--8<-- "examples/montecarlo.py:results"
```

```
                  count          mean           std           min            5%           50%           95%           max
GapAngle         2000.0  1.236226e+00      0.723843  6.096265e-02  2.983545e-01  1.087546e+00  2.675680e+00  3.140577e+00
GapHeight        2000.0  1.143269e-03      0.001238  9.070261e-08  4.384932e-05  7.442454e-04  3.596016e-03  1.128285e-02
BoltPreload      2000.0  2.732521e+06  27989.062233  2.600771e+06  2.688598e+06  2.732288e+06  2.777531e+06  2.916772e+06
Damage           2000.0  1.432734e-01      0.231241  1.589216e-03  1.232016e-02  7.530805e-02  4.771642e-01  4.232725e+00
AllowableDamage  2000.0  9.917610e-01      0.294008  2.884288e-01  5.983693e-01  9.423795e-01  1.534146e+00  2.654757e+00
FatigueLife      2000.0  7.505789e+02   1285.993299  8.399722e+00  5.808260e+01  3.838157e+02  2.471174e+03  2.216393e+04
Probability of failure over 30 years: 1.55e-02
```

!!! warning "Sample size"
    With 2000 realizations the estimate of a small probability is very
    uncertain: 31 failures out of 2000 give a coefficient of variation of
    about 18% on the estimate. For a target probability of failure of
    10<sup>-4</sup>, use at least 10<sup>6</sup> realizations, or a
    variance reduction technique.

This simulation does not apply a model uncertainty factor (stress
multiplication factor) to the stress ranges. Pass a constant or a sampler as
the `SMF` argument of `fatigue_case_sampler` to include it.
