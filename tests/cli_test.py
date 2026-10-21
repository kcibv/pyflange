import pytest
import os
import json
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import openpyxl
from typer.testing import CliRunner

from pyflange.cli import app
import pyflange.io as io

EXAMPLE_CASE_PATH = os.path.join(os.path.dirname(__file__), "..", "docsrc", "examples", "fatigue-case.json")

runner = CliRunner()


@pytest.fixture
def case_path (tmp_path):
    path = tmp_path / "case.json"
    shutil.copy(EXAMPLE_CASE_PATH, path)
    return str(path)


def invoke (*args):
    result = runner.invoke(app, [str(arg) for arg in args])
    assert result.exit_code == 0, result.output
    return result.output



def test_help ():
    output = invoke("--help")
    for command in ("determ", "sample", "inspect", "stats"):
        assert command in output
    assert "hist" not in output



class TestDeterm:

    def test_default (self, case_path, tmp_path):
        invoke("determ", case_path, "--odir", tmp_path)
        with open(tmp_path / "case-res.json") as file:
            res = json.load(file)
        fcase = io.load_fatiguecase(case_path)
        assert res['damage'] == pytest.approx(fcase.damage, rel=1e-12)
        assert res['FatigueLife'] == pytest.approx(fcase.fatigue_life, rel=1e-12)
        assert 'Point1' in res['FlangeSegment']

    def test_overrides (self, case_path, tmp_path):
        invoke("determ", case_path, "--odir", tmp_path, "--gap-angle", 45, "--gap-height", 1.5,
               "--tilt-angle", 0.2, "--smf", 1.1, "--xlsx")
        name = "case-res-GapA45.0deg-gaph1.5mm-TiltA0.2deg"
        with open(tmp_path / f"{name}.json") as file:
            res = json.load(file)
        assert res['FlangeSegment']['gap']['angle'] == pytest.approx(np.pi/4)
        assert res['FlangeSegment']['gap']['height'] == pytest.approx(0.0015)
        assert res['FlangeSegment']['tilt_angle'] == pytest.approx(0.2*np.pi/180)
        assert res['SMF'] == 1.1
        assert os.path.exists(tmp_path / f"{name}.xlsx")

    def test_custom_name (self, case_path, tmp_path):
        invoke("determ", case_path, "--odir", tmp_path, "--oname", "out.json")
        assert os.path.exists(tmp_path / "out.json")



class TestSampleInspectStats:

    def test_workflow (self, case_path, tmp_path):
        # Sample
        output = invoke("sample", case_path, "--size", 4, "--seed", 5, "--odir", tmp_path, "--xlsx")
        assert "seed 5" in output
        sample = pd.read_csv(tmp_path / "case-sample.csv", index_col=0)
        assert list(sample.columns) == io.SAMPLE_COLUMNS
        assert len(sample) == 4
        assert os.path.exists(tmp_path / "case-sample.xlsx")

        # Reproducibility
        invoke("sample", case_path, "--size", 4, "--seed", 5, "--odir", tmp_path, "--oname", "again.csv")
        pd.testing.assert_frame_equal(sample, pd.read_csv(tmp_path / "again.csv", index_col=0))

        # Inspect
        invoke("inspect", case_path, 2, "--odir", tmp_path, "--xlsx")
        with open(tmp_path / "case-realization-2.json") as file:
            res = json.load(file)
        assert res['damage'] == pytest.approx(sample.Damage[2], rel=1e-9)
        assert os.path.exists(tmp_path / "case-realization-2.xlsx")

        # Stats
        output = invoke("stats", tmp_path / "case-sample.csv")
        assert "MIN" in output
        assert "Damage" in output
        assert "sample size: 4" in output

    def test_inspect_custom_sample_path (self, case_path, tmp_path):
        invoke("sample", case_path, "--size", 2, "--seed", 5, "--odir", tmp_path, "--oname", "s.csv")
        invoke("inspect", case_path, 1, "--spath", tmp_path / "s.csv", "--odir", tmp_path, "--oname", "r.json")
        assert os.path.exists(tmp_path / "r.json")

    def test_stats_values (self, tmp_path):
        path = tmp_path / "sample.csv"
        pd.DataFrame({"X": [1.0, 2.0, 3.0]}).to_csv(path)
        output = invoke("stats", path)
        row = [line for line in output.splitlines() if "X" in line and "MIN" not in line][-1]
        values = [float(v) for v in row.split()[-7:]]
        assert values[0] == 1.0     # MIN (was mislabeled as MAX in Forrester)
        assert values[4] == 3.0     # MAX



def test_python_m_pyflange ():
    result = subprocess.run([sys.executable, "-m", "pyflange", "--help"], capture_output=True, encoding="utf-8",
                            env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    assert result.returncode == 0
    assert "determ" in result.stdout
