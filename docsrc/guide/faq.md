# FAQ and troubleshooting

## Installation

??? question "`'py' is not recognized as an internal or external command`"
    Python, or its *py launcher*, is not installed. See
    [Install PyFlange](../getting-started/install.md#3-check-that-python-works).

??? question "`'pyflange' is not recognized as an internal or external command`"
    The folder where pip installs commands is not in the `PATH`. Use
    `py -m pyflange` instead of `pyflange`: it always works.

??? question "`ModuleNotFoundError: No module named 'pyflange'`"
    PyFlange is not installed for the Python interpreter that runs your
    script. If you have several Python versions, install PyFlange with the
    same command you use to run the script, e.g. `py -3.12 -m pip install
    pyflange` for `py -3.12 script.py`. In VS Code, check the interpreter
    selected at the bottom-right of the window.

??? question "pip fails with a connection or SSL error"
    Your company network probably uses a proxy. Ask your IT department for
    its address and pass it to pip with `--proxy`, see
    [Install PyFlange](../getting-started/install.md#4-install-pyflange).

??? question "How do I know which version I have?"
    `py -m pip show pyflange`, or `py -m pyflange --version`. The changes
    between versions are listed in the [Changelog](../changelog.md).

## Modelling

??? question "`ValueError: Failure mode D detected, while only failure mode B is supported`"
    `validate()` found that the flange fails at ULS in a mode for which the
    polynomial model is not valid (see
    [Validity of the model](model.md#validity-of-the-model)). Typically the
    flange is too thin compared with the bolt (mode D or E) or the shell is
    too weak. Change the geometry, or use FEA.

??? question "Which shell stiffness should I use?"
    By default (`k_shell=None`) the simplified formula of ref. [2], eq. G.43
    is used, which is conservative. A better estimate is the FEA database
    of ref. [5], 9.3.1: pass `k_shell=shell_stiffness(R, s, gap_angle)`, or
    `"k_shell": "interp"` in a JSON case. The database covers diameters of 3
    to 10 m, shell thicknesses of 10 to 150 mm and gap angles of 10° to 180°;
    outside this range `shell_stiffness` returns `nan`.

??? question "Is `Fv` the preload before or after losses?"
    Before losses. The model subtracts the settlement and plastic losses
    defined by `DFT`, `settlement_factor` and `PPL`; the result is
    `fseg.bolt_preload`. If your preload already includes all the losses,
    leave `PPL` and `DFT` at zero. Note that the coefficient of variation
    term of ref. [2], eq. (23) is not applied by PyFlange: include it in
    `Fv`.

??? question "Why does the bolt force differ from my FEA results?"
    The polynomial model is an approximation calibrated on FEA results
    (refs. [1] and [5]), so some difference is expected; the stress
    multiplication factor of ref. [5], 6.5 accounts for this model
    uncertainty. First check that the gap, the
    shell stiffness, the inclination and the preload are the same in both
    models.

??? question "Can I model T-flanges in a JSON case file?"
    Not yet: the JSON format and the command line support L-flanges only. Use
    the Python API, see the [T-flange tutorial](../tutorials/tflange.md).

??? question "Should I use the AMD1 or the ED2 classes?"
    IEC 61400-6:2020+AMD1:2025 is the published standard: use the AMD1
    classes (`PolynomialLFlangeSegment`, `PolynomialTFlangeSegment`) for
    design. The ED2 classes implement a draft and are useful to anticipate
    the effect of the next edition.

## Results

??? question "The damage is zero or extremely small"
    The bolt force increases significantly only when the shell force
    approaches the gap opening region. Check that the Markov matrix contains
    **tower bending moments** in N·m (not kN·m or stresses), and that `Zg`
    is the dead weight **per bolt** (per bolt pair for T-flanges), in N.

??? question "The damage is `nan`"
    Usually an input is `nan`, for example `k_shell` from `shell_stiffness`
    outside the range of the FEA database. Check the intermediate quantities
    with `pyflange determ`, or print `fseg.shell_force_at_tensile_ULS` and
    the other points.

??? question "The Monte Carlo results change at every run"
    Use `--seed` on the command line, or `numpy.random.seed(...)` in Python.

## Getting help

If your question is not answered here, search the
[GitHub issues](https://github.com/kcibv/pyflange/issues) and open a new
one if needed. Include the PyFlange version, a minimal script or case file,
and the complete error message.
