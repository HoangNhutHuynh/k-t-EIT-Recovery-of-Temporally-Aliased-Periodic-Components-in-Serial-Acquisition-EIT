"""Figs. 3 to 5 of the k-t EIT manuscript: component images and error maps (IEEE two-column).
Data: *_data.mat written by run_export_fig_data.m (mesh, normalised nodal values, lung outlines, texts).

Part of k-t EIT (https://github.com/<user>/kt-eit). Author: Hoang Nhut Huynh. License: MIT.
Run from the repository root:  python figures/fig_images.py
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

import re, numpy as np, scipy.io as sio, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.tri import Triangulation

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5, "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.6})
COLS = ["Reference", "Naive", "Phase corr.", "FTC", "Temporal prior", "k-t (proposed)"]

def load(run, name):
    d = sio.loadmat(latest(run) + name + "_data.mat", squeeze_me=True, struct_as_record=False)
    P = np.atleast_1d(d["P"])
    rows = {}
    for p in P:          # group panels by row (axes y position), columns by x position
        rows.setdefault(round(float(p.pos[1]), 2), []).append(p)
    out = []
    for y in sorted(rows, reverse=True):
        out.append(sorted(rows[y], key=lambda p: p.pos[0]))
    return out, ListedColormap(d["cm"])

def texts(p):
    t = np.atleast_1d(p.txt)
    return [str(s) for s in t]

def metric(p, mode):
    for s in texts(p):
        if mode == "error":
            m = re.match(r"e=([\d.]+), L=([\d.]+)", s)
            if m: return float(m.group(1)), float(m.group(2))
        else:
            m = re.search(r"rho=([\d.]+)", s)
            if m: return float(m.group(1)), None
    return None

def draw_panel(ax, p, cmap, mode):
    V = np.asarray(p.V, float); F = np.asarray(p.F, int) - 1; C = np.asarray(p.C, float)
    tri = Triangulation(V[:, 0], V[:, 1], F)
    lo, hi = np.atleast_1d(p.clim).astype(float)
    ax.tripcolor(tri, C, shading="gouraud", cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
    lx, ly = np.atleast_1d(p.lx), np.atleast_1d(p.ly)
    lc = np.atleast_2d(p.lc)
    if lx.dtype != object:                        # a single line was squeezed to an array
        lx, ly = [lx], [ly]
    for x, y, c in zip(lx, ly, lc):
        x = np.asarray(x, float); y = np.asarray(y, float)
        grey = np.allclose(c, 0.5, atol=0.05)
        ax.plot(x, y, color="0.45" if grey else "k", lw=0.45 if grey else 0.7, solid_joinstyle="round")
    if np.size(p.mx):
        ax.plot(np.atleast_1d(p.mx), np.atleast_1d(p.my), marker="v", ms=3.2, color="k", ls="none")
    R = np.max(np.hypot(V[:, 0], V[:, 1])); lim = 1.08 * R
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off")

def figure(rows_spec, mode, fname, cbar_label):
    """rows_spec: list of (run, figname, row index in that figure, row label)."""
    nR, nC = len(rows_spec), 6
    W = 7.16; lm, rm = 0.46, 0.62; pw = (W - lm - rm) / nC
    ph = pw * 0.98; gap = 0.06; top = 0.20
    H = top + nR * (ph + gap) + 0.02
    fig = plt.figure(figsize=(W, H))
    cmap = None
    for r, (run, name, ri, label) in enumerate(rows_spec):
        rows, cmap = load(run, name)
        panels = rows[ri]
        assert len(panels) == nC, (run, name, len(panels))
        y0 = H - top - (r + 1) * ph - r * gap
        vals = [metric(p, mode) for p in panels]
        best = None
        cand = [(v[0], k) for k, v in enumerate(vals) if v is not None]
        if cand:
            best = (min(cand) if mode == "error" else max(cand))[1]
        for k, p in enumerate(panels):
            ax = fig.add_axes([(lm + k * pw) / W, y0 / H, pw / W, ph / H])
            draw_panel(ax, p, cmap, mode)
            if r == 0:
                ax.set_title(COLS[k], fontsize=7.5, fontweight="bold", pad=1.5)
            v = vals[k]
            if v is not None:
                s = (f"$e$={v[0]:.2f}, $e_L$={v[1]:.2f}" if mode == "error" else f"$\\rho$={v[0]:.2f}")
                ymin = np.min(np.asarray(p.V, float)[:, 1])
                ax.text(0.0, ymin - 0.03 * np.max(np.abs(np.asarray(p.V, float))), s, ha="center", va="top", fontsize=6.6,
                        fontweight="bold" if k == best else "normal")
        fig.text(0.06 / W, (y0 + ph / 2) / H, label, rotation=90, ha="left", va="center", fontsize=7, linespacing=1.15)
    # one colour bar for the whole figure
    cax = fig.add_axes([(W - rm + 0.10) / W, (H - top - nR * ph - (nR - 1) * gap) / H + 0.02,
                        0.07 / W, (nR * ph + (nR - 1) * gap) / H - 0.04])
    lo, hi = (0, 1) if mode == "error" else (-1, 1)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(lo, hi), cmap=cmap), cax=cax)
    cb.set_ticks([lo, (lo + hi) / 2, hi]); cb.set_ticklabels([f'{t:g}' for t in (lo, (lo + hi) / 2, hi)]); cb.ax.tick_params(labelsize=6.8, width=0.5, length=2)
    cb.outline.set_linewidth(0.5); cb.set_label(cbar_label, fontsize=7, labelpad=2)
    for ext in ("pdf", "svg"): fig.savefig(OUT / f"{fname}.{ext}", dpi=600)
    fig.savefig(OUT / f"{fname}.png", dpi=600)
    plt.close(fig)

def main():
    PIG = "*_E7decim_FINAL_E7_p1130122"
    HF, NEO, PL = "*_E7decim_FINAL_E7_hfov", "*_E7decim_FINAL_E7_neonate", "*_E7decim_E10_fig_g10-c005"

    # Fig. 3: PEEP-trial pig after surfactant, error maps of h5 to h8
    figure([(PIG, "fig_decim_error_main", i, l) for i, l in enumerate(
            ["$h_5$, 3.51 Hz", "$h_6$, 4.21 Hz", "$h_7$, 4.91 Hz", "$h_8$, 5.62 Hz"])],
           "error", "fig_pig_error", "Magnitude / max $|X_{\\mathrm{ref}}|$")
    CASES = [(HF, "HFOV pig\n$h_2$, 20.6 Hz"), (NEO, "Infant\ncardiac $h_2$, 4.9 Hz"), (PL, "Piglet 2\ncardiac, 3.5 Hz")]
    # Fig. 4: reconstructed images (real part after phase alignment)
    figure([(run, "fig_decim_aligned", 0, lab) for run, lab in CASES], "aligned", "fig_cases_images",
           "Re$(X)$ / max $|X_{\\mathrm{ref}}|$")
    # Fig. 5: error maps of the same cases
    figure([(run, "fig_decim_error", 0, lab) for run, lab in CASES], "error", "fig_cases_error",
           "Magnitude / max $|X_{\\mathrm{ref}}|$")


if __name__ == "__main__":
    main()
