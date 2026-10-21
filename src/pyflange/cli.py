# pyFlange - python library for large flanges design
# Copyright (C) 2024  KCI The Engineers B.V.,
#                     Siemens Gamesa Renewable Energy B.V.,
#                     Nederlandse Organisatie voor toegepast-natuurwetenschappelijk onderzoek TNO.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License version 3 for more details.
#
# You should have received a copy of the GNU General Public License
# version 3 along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""PyFlange command-line interface.

The command-line interface performs deterministic and probabilistic (Monte
Carlo) bolt fatigue analyses of flange cases described in the PyFlange JSON
format (see `pyflange.io`). It can be invoked as `pyflange` or as
`python -m pyflange`, followed by one of these commands:

- `determ`: deterministic fatigue analysis of a case.
- `sample`: Monte Carlo simulation of a case.
- `inspect`: rebuild and report a single Monte Carlo realization.
- `stats`: summary statistics of a Monte Carlo sample.

Run `pyflange --help` or `pyflange COMMAND --help` for the full list of options.
"""

import os
from typing import Annotated, Optional

from typer import Typer, Argument, Option, Exit


app = Typer(name="pyflange",
            help="Fatigue analysis of large bolted flanges, from PyFlange JSON case files.",
            no_args_is_help=True,
            add_completion=False)


def _version_callback(value: bool):
    if value:
        from importlib.metadata import version, PackageNotFoundError
        try:
            print(f"pyflange {version('pyflange')}")
        except PackageNotFoundError:
            print("pyflange (version unknown: package not installed)")
        raise Exit()


@app.callback()
def main(version: Annotated[Optional[bool], Option("--version", callback=_version_callback, is_eager=True,
                                                   help="Show the PyFlange version and exit.")] = None):
    """Fatigue analysis of large bolted flanges, from PyFlange JSON case files."""


def _output_path(odir, oname, default_name):
    return os.path.join(odir, oname or default_name)


def _write_dict(d, opath, xlsx):
    from .io import dict_to_json, dict_to_excel

    print(f"... writing data to json file '{opath}'")
    dict_to_json(d, opath)

    if xlsx:
        xlsx_opath = f"{os.path.splitext(opath)[0]}.xlsx"
        print(f"... writing data to excel file '{xlsx_opath}'")
        dict_to_excel(d, xlsx_opath)


@app.command()
def determ(
        cpath     : Annotated[str,   Argument(help="Path to the json file containing the case description.")],
        gap_angle : Annotated[Optional[float], Option(help="Custom gap angle in degrees. Defaults to input file value.")] = None,
        gap_height: Annotated[Optional[float], Option(help="Custom gap height in mm. Defaults to input file value.")] = None,
        tilt_angle: Annotated[Optional[float], Option(help="Custom tilt angle in degrees. Defaults to input file value.")] = None,
        smf       : Annotated[Optional[float], Option(help="Stress multiplication factor. Defaults to the input file SMF value, or 1.0.")] = None,
        odir      : Annotated[str,   Option(help="Path of the folder that will contain the output file(s). Defaults to cwd.")] = ".",
        oname     : Annotated[str,   Option(help="Custom name of the output json file, including extension. Optional.")] = "",
        xlsx      : Annotated[bool,  Option(help="If this flag is present, the output file will be saved in .xlsx format as well.")] = False
    ):
    """Performs a deterministic fatigue analysis.

    Given the path of a flange case input file, performs a deterministic
    fatigue analysis of the bolt, based on the data contained in the input
    file and writes the results to an output file.

    Some of the input parameters can be overridden directly from the
    command-line. When the gap angle is overridden, the shell stiffness is
    re-interpolated for the new angle.
    """
    from math import pi
    from .io import load_inputfile, upgrade_case, create_fatiguecase, fatiguecase_to_dict
    from .flangesegments import shell_stiffness

    deg = pi/180
    mm = 0.001

    print(f"Performing deterministic fatigue analysis for case '{cpath}' ...")
    print(f"... loading data from '{cpath}'")
    params = upgrade_case(load_inputfile(cpath))

    print(f"... building the model")
    fseg_params = params['FlangeSegment']
    suffix = ""

    if gap_angle is not None:
        fseg_params['gap']['angle'] = gap_angle*deg
        fseg_params['k_shell'] = float(shell_stiffness(fseg_params['R'], fseg_params['s'], fseg_params['gap']['angle']))
        suffix += f"-GapA{gap_angle}deg"

    if gap_height is not None:
        fseg_params['gap']['height'] = gap_height*mm
        suffix += f"-gaph{gap_height}mm"

    if tilt_angle is not None:
        fseg_params['tilt_angle'] = tilt_angle*deg
        suffix += f"-TiltA{tilt_angle}deg"

    fcase = create_fatiguecase(params, SMF=smf)

    opath = _output_path(odir, oname, f"{os.path.splitext(os.path.basename(cpath))[0]}-res{suffix}.json")
    _write_dict(fatiguecase_to_dict(fcase, verbose=True), opath, xlsx)

    print(f"... DONE!")


@app.command()
def sample(
        fpath: Annotated[str,  Argument(help="Path of the file containing the fatigue case data.")],
        size : Annotated[int,  Option(help="Number of realizations to be generated.")] = 100000,
        seed : Annotated[Optional[int], Option(help="Seed of the random number generator, for reproducible simulations. Optional.")] = None,
        odir : Annotated[str,  Option(help="Path of the folder that will contain the csv file to which the sample will be written.")] = ".",
        oname: Annotated[str,  Option(help="Custom name of the output csv file, including extension. Optional.")] = "",
        xlsx : Annotated[bool, Option(help="If this flag is present, the output file will be saved in .xlsx format as well.")] = False
    ):
    """Runs a Monte Carlo simulation and saves the results to a file.

    The output CSV file contains one row per realization, with the random
    input parameters and the resulting damage and fatigue life. Values are
    in SI units (m, N, rad, Pa), fatigue life in the units of the Markov
    matrix duration.
    """
    from .io import load_inputfile, sample_fatiguecase

    full_path = os.path.abspath(fpath)

    print(f"Performing probabilistic analysis for case {fpath} ...")
    print(f"... loading data from '{full_path}'")
    data = load_inputfile(full_path)

    seed_info = f" with seed {seed}" if seed is not None else ""
    print(f"... generating random sample of {size} realizations{seed_info}")

    def progress(i):
        if i % 1000 == 0:
            print(f"... ... generated {i} realizations")

    sampdf = sample_fatiguecase(data, size, seed=seed, progress=progress)

    opath = _output_path(odir, oname, f"{os.path.splitext(os.path.basename(fpath))[0]}-sample.csv")
    print(f"... writing results to file '{opath}'")
    sampdf.to_csv(opath)

    if xlsx:
        xlsx_opath = f"{os.path.splitext(opath)[0]}.xlsx"
        print(f"... writing results to file '{xlsx_opath}'")
        sampdf.to_excel(xlsx_opath, sheet_name="Sample")

    print(f"... DONE!")


@app.command()
def inspect(
        cpath: Annotated[str,  Argument(help="Path to the json file containing the case description.")],
        real : Annotated[int,  Argument(help="Number of the realization to be inspected.")],
        odir : Annotated[str,  Option(help="Path of the folder that will contain the json file to which the data will be written.")] = ".",
        oname: Annotated[str,  Option(help="Custom name of the output json file, including extension. Optional.")] = "",
        spath: Annotated[str,  Option(help="Path to the CSV file containing the sample. Defaults to '<case>-sample.csv'.")] = "",
        xlsx : Annotated[bool, Option(help="If this flag is present, the output file will be saved in .xlsx format as well.")] = False
    ):
    """Extracts a Monte Carlo realization for inspection.

    Rebuilds the fatigue case of a given realization of a sample generated
    with the `sample` command, and writes its detailed report to a file.
    """
    import pandas as pd
    from .io import load_inputfile, create_realization, fatiguecase_to_dict

    print(f"Inspecting realization number {real} from case {cpath} ...")
    print(f"... loading data from '{cpath}'")
    data = load_inputfile(cpath)

    spath = spath or f"{os.path.splitext(cpath)[0]}-sample.csv"
    print(f"... loading sample from '{spath}'")
    sampdf = pd.read_csv(spath, index_col=0)

    print(f"... re-assembling data of realization number {real}")
    fcase = create_realization(data, sampdf.loc[real])

    opath = _output_path(odir, oname, f"{os.path.splitext(os.path.basename(cpath))[0]}-realization-{real}.json")
    _write_dict(fatiguecase_to_dict(fcase, verbose=True), opath, xlsx)

    print(f"... DONE!")


@app.command()
def stats(fpath: Annotated[str, Argument(help="Path of the csv file containing the sample data.")]):
    """Prints statistics of each parameter of a sample.

    Given a sample file path, prints synthetic statistics (minimum, 5%
    quantile, average, 95% quantile, maximum, standard deviation and
    coefficient of variation) of all the parameters in the sample.
    """
    import numpy as np
    import pandas as pd

    print(f"Generate statistics for sample in '{fpath}' ...")

    print(f"... loading sample from '{fpath}'")
    sampdf = pd.read_csv(fpath, index_col=0)

    print(f"... computing statistics")
    stats_df = pd.DataFrame({
        f"... ... {name}": [
            np.min(series),
            np.quantile(series, 0.05),
            np.average(series),
            np.quantile(series, 0.95),
            np.max(series),
            np.std(series),
            np.std(series) / np.average(series)
        ]
        for name, series in sampdf.items()
    }, index=['MIN', 'P05', 'AVG', 'P95', 'MAX', 'STD', 'CoV'])

    print(f"... showing results\n")
    print(stats_df.transpose().to_string())
    print(f"\n... ... sample size: {len(sampdf)}")
    print(  f"... ... units of measurement: m, N, rad, Pa, fatigue life in units of the Markov matrix duration.")

    print(f"\n... DONE!")


if __name__ == "__main__":
    app()
