"""PyFlange: a Python library for the design of large bolted flange connections.

The library is organized in the following modules:

- `pyflange.bolts`: bolts, washers and nuts.
- `pyflange.flangesegments`: mechanical models of single-bolt flange segments
  (L-flanges and T-flanges).
- `pyflange.fatigue`: Markov matrices, S-N curves and bolt fatigue analysis.
- `pyflange.stats`: probability distributions and random samplers for
  probabilistic (Monte Carlo) analyses.
- `pyflange.io`: reading and writing fatigue cases in the PyFlange JSON format.
- `pyflange.cli`: command-line interface (`pyflange` or `python -m pyflange`).
- `pyflange.utils`: internal helpers.
- `pyflange.gap`: legacy module, kept for backwards compatibility.
"""
