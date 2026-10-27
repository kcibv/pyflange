# Getting started

These pages take you from a Windows computer without Python to your first
flange fatigue analysis. No programming experience is needed: every command
is given in full, and every step tells you what you should see.

There are two ways of using PyFlange:

<div class="grid cards" markdown>

-   :material-console:{ .lg .middle } **Without writing code**

    ---

    Describe the flange in a text file (a *case file*) and run the analysis
    with one command. Results are written to JSON, CSV or Excel files.

    Good for routine checks and for Monte Carlo simulations of a single
    flange.

-   :material-language-python:{ .lg .middle } **With Python scripts**

    ---

    Build the flange model in a short Python script. You have access to every
    intermediate quantity and you can run parameter studies, plots and
    batches.

    Good for design iterations and for anything the command line does not
    cover (for example T-flanges).

</div>

Both ways use the same model and give the same results.

## The steps

1. [**Install PyFlange**](install.md): install Python, open a terminal and
   install the package. About 10 minutes.
2. [**First analysis without code**](first-analysis.md): download an example
   case file, run a fatigue analysis and open the results in Excel.
3. [**First Python script**](first-script.md): write, run and modify a short
   script that computes a bolt force curve.

When you are done, continue with the [tutorials](../tutorials/index.md).

!!! tip "Already familiar with Python?"
    Install the package with `pip install pyflange` and go straight to the
    [tutorials](../tutorials/index.md).
