# Install PyFlange

This page explains how to install Python and PyFlange on Windows 10 or 11.
You do not need administrator rights.

## 1. Install Python

PyFlange needs **Python 3.9 or later**. If you are not sure whether Python is
already installed, skip to [step 3](#3-check-that-python-works): if the check
works, you can skip the installation.

1. Go to <https://www.python.org/downloads/> and click the yellow
   **Download Python 3.x.x** button.
2. Run the downloaded file (`python-3.x.x-amd64.exe`).
3. In the first window of the installer:
    - tick **Add python.exe to PATH** at the bottom of the window;
    - click **Install Now**.
4. Wait until the installer shows *Setup was successful*, then click
   **Close**.

!!! note "Company computers"
    If your company installs software for you (e.g. through a software
    center), ask for *Python 3.12* or later and check that the **`py`**
    command works as described in step 3.

## 2. Open a terminal

A *terminal* (also called *command prompt*) is a window where you type
commands instead of clicking.

1. Press the ++win++ key, type `cmd` and press ++enter++.
2. A black window opens, showing a line like:

    ```
    C:\Users\yourname>
    ```

    This is the *prompt*. It shows the current folder, here your user
    folder. You type a command after the `>` and press ++enter++ to run it.

!!! tip "Windows Terminal and PowerShell"
    You can also use *Windows Terminal* or *PowerShell*. All the commands in
    this documentation work in all of them.

A few commands you will use:

| Command | What it does |
|---------|--------------|
| `cd Documents` | Go to the `Documents` folder inside the current folder. |
| `cd ..` | Go up one folder. |
| `cd C:\Projects\flange` | Go to a specific folder. Use quotes if the path contains spaces: `cd "C:\My Projects"`. |
| `mkdir flange` | Create a folder called `flange`. |
| `dir` | List the files in the current folder. |
| `start .` | Open the current folder in the File Explorer. |
| `cls` | Clear the window. |

Press ++up++ to recall the previous command, and ++tab++ to complete a
file or folder name.

## 3. Check that Python works

In the terminal, type the following command and press ++enter++:

```doscon
C:\Users\yourname>py --version
Python 3.12.7
```

You should see `Python` followed by a version number of 3.9 or higher.

??? failure "`'py' is not recognized as an internal or external command`"
    Python is not installed, or it was installed without the *py launcher*.
    Run the installer again: choose **Modify**, make sure **py launcher** is
    selected and complete the installation. Then close the terminal and
    open a new one: a terminal only sees programs that were installed
    before it was opened.

??? failure "The Microsoft Store opens"
    Windows has a shortcut that opens the Store when you type `python`.
    Use `py` instead of `python`, as in this documentation, or install
    Python from python.org as described above.

## 4. Install PyFlange

Type:

```doscon
C:\Users\yourname>py -m pip install pyflange
```

`pip` is the Python package installer. It downloads PyFlange and the
packages it needs (`numpy`, `scipy`, `pandas`, ...). After a few seconds to a
few minutes, the last line should be:

```
Successfully installed ... pyflange-x.y.z ...
```

!!! note "Behind a company proxy"
    If pip cannot connect to the internet, ask your IT department for the
    proxy address and use:
    `py -m pip install --proxy http://proxy.company.com:8080 pyflange`

## 5. Check that PyFlange works

```doscon
C:\Users\yourname>py -m pyflange --version
pyflange x.y.z
```

You can now use the PyFlange command-line tool. You can always start it as
`py -m pyflange`. Usually the shorter form `pyflange` works too:

```doscon
C:\Users\yourname>pyflange --help
```

??? failure "`'pyflange' is not recognized ...`"
    The folder where pip installs commands is not in your `PATH`. pip
    prints a warning about it during the installation. You don't need to
    fix it: just use `py -m pyflange` instead of `pyflange`.

## Updating PyFlange

To install the latest version:

```doscon
C:\Users\yourname>py -m pip install --upgrade pyflange
```

To check which version you have: `py -m pip show pyflange`.

??? info "Optional: a separate environment per project"
    Engineers who use Python for several projects often create a *virtual
    environment* per project, so that each project keeps its own package
    versions. In the project folder:

    ```doscon
    C:\Projects\flange>py -m venv .venv
    C:\Projects\flange>.venv\Scripts\activate
    (.venv) C:\Projects\flange>pip install pyflange
    ```

    Run `.venv\Scripts\activate` every time you open a new terminal for the
    project. This is not required to use PyFlange.

## Next step

[Run your first analysis :octicons-arrow-right-24:](first-analysis.md)
