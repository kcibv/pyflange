# Tests that keep the documentation in sync with the code:
# - the tutorial scripts in docsrc/examples and the README examples run without errors;
# - every command and option of the command-line interface is documented.

import os
import re
import shutil
import runpy

import pytest

DOCSRC = os.path.join(os.path.dirname(__file__), "..", "docsrc")
EXAMPLES = os.path.join(DOCSRC, "examples")


def example_scripts():
    return sorted(name for name in os.listdir(EXAMPLES) if name.endswith(".py"))


@pytest.mark.parametrize("script", example_scripts())
def test_example_runs(script, tmp_path, monkeypatch, capsys):
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        pass

    # Run in an empty folder, with the example case file available
    shutil.copy(os.path.join(EXAMPLES, "fatigue-case.json"), tmp_path)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PYFLANGE_SAMPLE_SIZE", "20")

    runpy.run_path(os.path.join(EXAMPLES, script), run_name="__main__")


def test_readme_examples_run():
    readme_path = os.path.join(os.path.dirname(__file__), "..", "README.md")
    with open(readme_path, encoding="utf-8") as file:
        blocks = re.findall(r"```python\n(.*?)```", file.read(), re.S)
    assert blocks
    namespace = {}
    for block in blocks:        # the blocks build on each other
        exec(block, namespace)


def test_cli_documented():
    import typer.main
    from pyflange.cli import app

    with open(os.path.join(DOCSRC, "cli.md"), encoding="utf-8") as file:
        cli_doc = file.read()

    command = typer.main.get_command(app)
    for name, subcommand in command.commands.items():
        assert f"### {name}" in cli_doc, f"Command '{name}' is not documented in docsrc/cli.md"
        for param in subcommand.params:
            for opt in param.opts:
                if opt.startswith("--") and opt != "--help":
                    assert re.search(rf"`{opt}`", cli_doc), \
                        f"Option '{opt}' of '{name}' is not documented in docsrc/cli.md"
