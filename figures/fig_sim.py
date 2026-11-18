"""Fig. 2 of the k-t EIT manuscript: principal angles and simulation sweeps (IEEE two-column width).

Part of k-t EIT (https://github.com/<user>/kt-eit). Author: Hoang Nhut Huynh. License: MIT.
Run from the repository root:  python figures/fig_sim.py
Reads the MATLAB run logs in runs/ (newest run of each kind); writes the figure next to this script.
"""
import sys, glob, os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"
OUT = Path(__file__).resolve().parent

def latest(pattern):
    """Newest run folder in runs/ whose name matches pattern (e.g. '*_E3collision_FINAL_E3')."""
    c = sorted(glob.glob(str(RUNS / pattern)))
    if not c:
        sys.exit(f"no run matching {pattern} in {RUNS}: run the MATLAB scripts first (see README)")
    return c[-1] + os.sep

import csv, collections, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def load(pattern): return list(csv.DictReader(open(latest(pattern) + "metrics.csv")))

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "lines.linewidth": 1.3, "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False})
STY = {  # one colour per method, shared by all figures of the paper
    "naive":     ("Naive",              "#7a7a7a", "s", (0, (1, 1.2))),
    "yerworth":  ("Phase corr.",        "#2f6db5", "^", (0, (4, 1.5))),
    "ftc":       ("FTC",                "#e08a1e", "v", (0, (2, 1, 0.6, 1))),
    "kt_single": ("k-t, single",        "#8c5bb5", "D", "-"),
    "kt":        ("k-t, joint (proposed)", "#c62f2f", "o", "-"),
}
fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.5), gridspec_kw={"wspace": 0.42})

# (a) principal angles
A = collections.defaultdict(dict)
for x in load("*_E3collision_FINAL_E3"):
    if x["metric"] == "min_angle_deg" and x["param"] == "sep_per_stim_duty":
        A[float(x["pvalue"])][int(x["method"][4:])] = float(x["value"])
cols = plt.cm.viridis([0.08, 0.32, 0.55, 0.74, 0.90])
for c, d in zip(cols, sorted(A)):
    k = sorted(A[d]); ax[0].plot(k, [A[d][i] for i in k], color=c, lw=1.3, label=f"{d:g}")
ax[0].axvline(104, color="0.6", lw=0.7, ls=":")
ax[0].text(108, 64, "$M/2$", fontsize=7, color="0.35", ha="left", va="top")
ax[0].legend(title="Duty $\\delta$", fontsize=7, title_fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.62), handlelength=1.4, labelspacing=0.25)
ax[0].set_xlim(0, 200); ax[0].set_ylim(0, 68)
ax[0].set_xlabel("Subspace rank $r$"); ax[0].set_ylabel("Min. principal angle (deg)")

import numpy as np
NL = collections.defaultdict(list)
for x in csv.DictReader(open(RUNS / "R1_X2_nonlinear_sim.csv")):
    if x["metric"] == "rho_true" and x["component"] == "hfov_h3":
        NL[(x["sweep"], float(x["value"]), float(x["snr_db"]), x["method"])].append(float(x["metric_value"]))
SNR = 30.0
def draw(a, sweep, methods):
    xs = sorted({k[1] for k in NL if k[0] == sweep})
    for m in methods + ["simref"]:
        mu = np.array([np.mean(NL[(sweep, v, SNR, m)]) for v in xs]); sd = np.array([np.std(NL[(sweep, v, SNR, m)]) for v in xs])
        if m == "simref":
            a.plot(xs, mu, color="k", lw=0.8, ls="--", zorder=2); continue
        lab, c, mk, ls = STY[m]
        a.fill_between(xs, mu - sd, mu + sd, color=c, alpha=0.18, lw=0, zorder=1)
        a.plot(xs, mu, color=c, marker=mk, ms=3.2, ls=ls, lw=1.6 if m == "kt" else 1.1,
               zorder=4 if m == "kt" else 3, mec="white", mew=0.4)
    if sweep == "sim.f0":   # k-t joint at 20 dB for reference
        mu20 = [np.mean(NL[(sweep, v, 20.0, "kt")]) for v in xs]
        a.plot(xs, mu20, color=STY["kt"][1], lw=0.8, ls=(0, (1, 1)), zorder=3)

# (b) wrong assumed duty
draw(ax[1], "rec.timing_assumed.duty", ["naive", "yerworth", "ftc", "kt"])
ax[1].axvline(0.9, color="0.6", lw=0.7, ls=":"); ax[1].text(0.885, 0.06, "true", fontsize=7, color="0.35", ha="right")
ax[1].set_xlabel("Assumed duty $\\hat{\\delta}$"); ax[1].set_ylabel("$\\rho_{\\mathrm{true}}$, 3rd harmonic")
ax[1].set_ylim(0, 1.02); ax[1].set_xlim(0.48, 1.02)

# (c) sweep of f0
ax[2].axvspan(7.33, 20.5, color="0.94", zorder=0)
draw(ax[2], "sim.f0", ["naive", "yerworth", "ftc", "kt_single", "kt"])
ax[2].annotate("$3f_0$ folds onto $f_0$", xy=(11.25, 0.62), xytext=(13.0, 0.62), va="center", fontsize=6.8, color="0.25",
               arrowprops=dict(arrowstyle="-", lw=0.6, color="0.4"))
ax[2].text(20.3, 0.80, "k-t, 20 dB", fontsize=6.3, color=STY["kt"][1], ha="right", va="top")
ax[2].set_xlabel("Fundamental $f_0$ (Hz)"); ax[2].set_ylabel("$\\rho_{\\mathrm{true}}$, 3rd harmonic")
ax[2].set_ylim(0, 1.02); ax[2].set_xlim(5.5, 20.5)

for a, l in zip(ax, "abc"):
    a.text(-0.30, 1.04, f"({l})", transform=a.transAxes, fontsize=8.5, fontweight="bold", va="bottom")
h = [Line2D([], [], color=STY[m][1], marker=STY[m][2], ms=3.2, ls=STY[m][3], lw=1.3, mec="white", mew=0.4,
            label=STY[m][0]) for m in STY] + [Line2D([], [], color="k", ls="--", lw=0.8, label="Simultaneous acq.")]
fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.66, 1.0), ncol=3, fontsize=7.2,
           columnspacing=1.2, handlelength=2.4, handletextpad=0.5)
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.18, top=0.76)
for ext in ("pdf", "svg"): fig.savefig(OUT / ("fig_sim." + ext))
fig.savefig(OUT / "fig_sim.png", dpi=600)
