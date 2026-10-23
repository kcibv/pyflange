# Generates the figures of the documentation, in a light and a dark version.
#
# Run from the repository root after changes to the model or to the
# tutorial examples:
#
#     py docsrc/assets/img/make_figures.py

import os
import runpy
import contextlib, io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "..", "..", "examples")

THEMES = {
    "light": dict(fg="#1f2933", grid="#d5dbe3", c1="#1f63b0", c2="#c99700"),
    "dark":  dict(fg="#d6dde7", grid="#2e3949", c1="#79b0ef", c2="#ffd43b"),
}


def run_example(name):
    # Runs a tutorial script in a temporary folder and returns its globals
    cwd = os.getcwd()
    os.chdir(os.environ.get("TEMP", "/tmp"))
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return runpy.run_path(os.path.join(EXAMPLES, name))
    finally:
        os.chdir(cwd)


def style(ax, th):
    ax.set_facecolor("none")
    ax.tick_params(colors=th["fg"])
    for spine in ax.spines.values():
        spine.set_color(th["grid"])
    ax.xaxis.label.set_color(th["fg"])
    ax.yaxis.label.set_color(th["fg"])
    ax.title.set_color(th["fg"])
    ax.grid(True, color=th["grid"], linewidth=0.6)


def transfer_functions():
    import numpy as np
    g = run_example("lflange.py")
    fseg = g["fseg"]
    kN = 1000
    Z = np.linspace(-700*kN, 2000*kN, 400)
    points = [
        ("P1", fseg.shell_force_at_rest),
        ("P2", fseg.shell_force_at_tensile_ULS),
        ("P3", fseg.shell_force_at_small_displacement),
        ("P4", fseg.shell_force_at_closed_gap),
    ]
    for name, th in THEMES.items():
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))
        fig.patch.set_alpha(0)
        ax1.plot(Z/kN, fseg.bolt_axial_force(Z)/kN, color=th["c1"], lw=2)
        ax2.plot(Z/kN, fseg.bolt_bending_moment(Z)/1000, color=th["c1"], lw=2)
        for label, Zp in points:
            Fs = fseg.bolt_axial_force(Zp)/kN
            Ms = fseg.bolt_bending_moment(Zp)/1000
            ax1.plot(Zp/kN, Fs, "o", color=th["c2"])
            ax1.annotate(label, (Zp/kN, Fs), textcoords="offset points", xytext=(-6, 8),
                         color=th["fg"], fontsize=9)
            ax2.plot(Zp/kN, Ms, "o", color=th["c2"])
            ax2.annotate(label, (Zp/kN, Ms), textcoords="offset points", xytext=(-6, 8),
                         color=th["fg"], fontsize=9)
        ax1.set_xlabel("Shell force Z [kN]")
        ax1.set_ylabel("Bolt axial force Fs [kN]")
        ax2.set_xlabel("Shell force Z [kN]")
        ax2.set_ylabel("Bolt bending moment Ms [kNm]")
        style(ax1, th)
        style(ax2, th)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, f"lflange-transfer-functions-{name}.png"), dpi=130, transparent=True)
        plt.close(fig)


if __name__ == "__main__":
    transfer_functions()
