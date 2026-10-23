# Writing case files from Python

The [JSON case files](../json-format.md) used by the command line can also be
read, modified and written from Python with [`pyflange.io`](../apis/io.md).
This is convenient for running variants of a case and for collecting the
results of many cases in one table.

**Download the full script:**
[case_file.py](../examples/case_file.py){ download="case_file.py" }. Run it
in a folder that contains
[fatigue-case.json](../examples/fatigue-case.json){ download="fatigue-case.json" }.

## 1. Load a case

```python
--8<-- "examples/case_file.py:load"
```

```
Damage of the example case: 0.501
```

[`load_inputfile`][pyflange.io.load_inputfile] returns the content of the file
as a dictionary, with the same structure as the JSON file.
[`create_fatiguecase`][pyflange.io.create_fatiguecase] builds the
[`BoltFatigueAnalysis`][pyflange.fatigue.BoltFatigueAnalysis] from it. To do
both in one step, use [`load_fatiguecase`][pyflange.io.load_fatiguecase].

## 2. Write variants of the case

Change the dictionary and save it with
[`dict_to_json`][pyflange.io.dict_to_json]. Here we create one case per
design gap angle, with the AMD1 design gap height of each angle:

```python
--8<-- "examples/case_file.py:variants"
```

Since the example case has `"k_shell": "interp"`, the shell stiffness is
interpolated for the new gap angle when the case is loaded.

## 3. Export a case with its results

[`fatiguecase_to_dict`][pyflange.io.fatiguecase_to_dict] converts an analysis
back to the JSON format. With `verbose=True` it adds the results and the
intermediate quantities of the model, like the output of
`pyflange determ`:

```python
--8<-- "examples/case_file.py:export"
```

The same function works for an analysis built in Python, for example the
`fatigue` object of the [L-flange tutorial](lflange.md): this is the way to
create a case file from a Python model.

```python
from pyflange.io import fatiguecase_to_dict, dict_to_json
dict_to_json(fatiguecase_to_dict(fatigue), "my-case.json")
```

!!! note
    The JSON format currently supports L-flange segments only.

## 4. Run a batch

```python
--8<-- "examples/case_file.py:batch"
```

```
            case  gap angle [deg]   damage  fatigue life [years]
case-gap120.json            120.0 1.080555             27.763512
 case-gap30.json             30.0 0.503086             59.632003
 case-gap60.json             60.0 0.763320             39.302010
 case-gap90.json             90.0 1.169948             25.642169
```

The table is also saved to `results.xlsx`. To add a Monte Carlo simulation of
each case, use [`sample_fatiguecase`][pyflange.io.sample_fatiguecase], which
returns the same table as the `sample` command.
