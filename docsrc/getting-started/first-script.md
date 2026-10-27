# First Python script

In this step you write a short Python script that builds an L-flange segment
and prints its bolt force curve. You need PyFlange [installed](install.md).

## 1. Open IDLE

Python comes with a simple editor called **IDLE**. Press ++win++, type
`idle` and press ++enter++. A window titled *IDLE Shell* opens.

!!! tip "Other editors"
    For larger projects, [Visual Studio Code](https://code.visualstudio.com/)
    with the Python extension is a better editor. Everything on this page
    works the same way there: write the script, save it, run it.

## 2. Write the script

In IDLE, choose **File > New File**. An empty editor window opens. Copy the
following script into it (use the copy button at the top-right of the box):

```python title="first_script.py"
--8<-- "examples/first_script.py"
```

Choose **File > Save As...**, go to your `flange-analyses` folder and save
the file as `first_script.py`.

## 3. Run it

Press ++f5++ (or **Run > Run Module**). The output appears in the IDLE Shell
window:

```
Shell force   -200 kN  ->  bolt force  2876.0 kN
Shell force      0 kN  ->  bolt force  2878.2 kN
Shell force    200 kN  ->  bolt force  2892.4 kN
Shell force    400 kN  ->  bolt force  2919.9 kN
Shell force    600 kN  ->  bolt force  2960.7 kN
Shell force    800 kN  ->  bolt force  3014.9 kN
Shell force   1000 kN  ->  bolt force  3082.3 kN
```

You can also run the script from the terminal:

```doscon
C:\Users\yourname\flange-analyses>py first_script.py
```

## 4. Understand the script

- **`from ... import ...`** lines load the PyFlange objects used by the
  script: a standard bolt and nut, the L-flange segment model and the gap.
- **Units.** PyFlange does not know about units: it expects all values in one
  consistent system. The script uses meters, newtons and radians, and defines
  `mm`, `kN` and `deg` so that you can write `150*mm` instead of `0.150`.
  See [Units and conventions](../guide/units.md).
- **`PolynomialLFlangeSegment(...)`** creates the model of one flange
  *segment*: the slice of the flange that contains one bolt. The symbols `a`,
  `b`, `s`, `t` follow IEC 61400-6 and are explained in the
  [Nomenclature](../guide/nomenclature.md). `Zg` is the shell force per bolt
  due to the dead weight (negative: compression).
- **`fseg.bolt_axial_force(Z)`** returns the bolt force for a shell force
  `Z` acting on the segment.
- **`print(f"...")`** prints a line of text. The expressions between braces
  are replaced by their value; `:7.1f` means *7 characters wide, 1 decimal*.

## 5. Change something

Try these changes one at a time. Save (++ctrl+s++) and run (++f5++) after
each one.

- Make the gap larger: `gap = Gap(height=2.0*mm, angle=30*deg)`. The bolt
  force grows faster with the shell force, because a larger part of the
  load goes through the bolt until the gap closes.
- Reduce the preload: `Fv = 2500*kN`.
- Print the bolt bending moment as well: add the line
  `print(fseg.bolt_bending_moment(Z))` inside the `for` loop (indented like
  the other lines of the loop).

??? failure "`IndentationError` or `SyntaxError`"
    Python uses indentation to group lines. The lines inside the `for` loop
    must be indented by the same number of spaces (4 in this script). A
    missing comma or parenthesis also causes a `SyntaxError`: the error
    message shows the line where Python got confused.

??? failure "`ModuleNotFoundError: No module named 'pyflange'`"
    PyFlange is not installed for the Python you are using. Run
    `py -m pip install pyflange` again, and check that IDLE was started
    from the same Python installation (`py -m idlelib` starts it).

## Next steps

You now know how to build a flange segment in Python. The
[L-flange tutorial](../tutorials/lflange.md) continues from here: it checks
the validity of the model, plots the bolt force and moment curves and
performs the fatigue analysis.
