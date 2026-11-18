"""Graphical abstract of the k-t EIT manuscript (aspect 660:295). Numbers in panel 4 are those of Table IV.

Part of k-t EIT (https://github.com/<user>/kt-eit). Author: Hoang Nhut Huynh. License: MIT.
Run from the repository root:  python figures/graphical_abstract.py
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

import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, FancyBboxPatch, Polygon, Rectangle, FancyArrowPatch, Arc
sys.path.insert(0, str(OUT))
import fig_images as FI   # reuse the exported mesh data of Fig. 4

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "svg.fonttype": "none", "mathtext.fontset": "dejavusans"})
INK, GREY, LIGHT, RED, BLUE = "#1f2937", "#6b7280", "#eef1f5", "#c62f2f", "#2f6db5"
PIG, PIGL, SKIN = "#f4b6a6", "#e79c88", "#f5d1ad"
W, H = 7.33, 3.277            # 660 x 295 proportion
fig = plt.figure(figsize=(W, H))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis("off")

# ---------- frame of the four steps
ws = [1.62, 1.66, 1.62, 2.05]; gap = 0.09; xs = [0.08]
for w_ in ws[:-1]: xs.append(xs[-1] + w_ + gap)
y0, hh = 0.42, 2.62
titles = ["Problem: aliasing", "k-t forward model", "Real data", "Result"]
for k, (x, w, t) in enumerate(zip(xs, ws, titles)):
    bg.add_patch(FancyBboxPatch((x, y0), w, hh, boxstyle="round,pad=0,rounding_size=0.08", fc="white", ec="#c9cfd8", lw=0.8))
    bg.add_patch(Circle((x + 0.15, y0 + hh - 0.17), 0.09, fc=INK, ec="none"))
    bg.text(x + 0.15, y0 + hh - 0.17, str(k + 1), color="white", ha="center", va="center", fontsize=7.5, fontweight="bold")
    bg.text(x + 0.30, y0 + hh - 0.17, t, ha="left", va="center", fontsize=7.6, fontweight="bold", color=INK)
    if k < 3:
        bg.add_patch(FancyArrowPatch((x + w + 0.01, y0 + hh / 2), (xs[k + 1] - 0.01, y0 + hh / 2),
                                     arrowstyle="-|>", mutation_scale=10, color="#9aa3b2", lw=1.2))
bg.text(W / 2, 0.20, "When the measurements of a frame are acquired serially, the frame rate no longer limits a periodic component of known frequency.",
        ha="center", va="center", fontsize=7.0, color=INK, style="italic")

# ---------- 1: problem (samples of one measurement, true vs folded frequency)
a1 = fig.add_axes([(xs[0] + 0.08) / W, (y0 + 0.30) / H, (ws[0] - 0.16) / W, 1.85 / H]); a1.axis("off")
a1.set_xlim(-0.05, 3.05); a1.set_ylim(-1.75, 2.35)
t = np.linspace(0, 3, 600); tau = 0.40
a1.plot(t, 0.95 + 0.55 * np.cos(2 * np.pi * 1.25 * t), color=RED, lw=1.3)
a1.plot(t, 0.95 + 0.55 * np.cos(2 * np.pi * 0.25 * t + 2 * np.pi * tau), color=GREY, lw=1.1, ls=(0, (3, 1.6)))
ts = np.arange(3) + tau
a1.plot(ts, 0.95 + 0.55 * np.cos(2 * np.pi * 1.25 * ts), "o", color=INK, ms=3.6, zorder=5)
for f in range(3):
    a1.add_patch(Rectangle((f, -0.62), 1, 0.52, fc=LIGHT, ec="#b6bfcc", lw=0.6))
    for s in range(8):
        a1.add_patch(Rectangle((f + 0.03 + s * 0.115, -0.56), 0.085, 0.40, fc="#1f4f8f" if s == 3 else "#8fb2dd", ec="none"))
    a1.plot([f + tau, f + tau], [-0.10, 0.95 + 0.55 * np.cos(2 * np.pi * 1.25 * (f + tau))], color="#b6bfcc", lw=0.6, ls=":")
a1.annotate("", xy=(tau, -0.80), xytext=(0, -0.80), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=INK, mutation_scale=6))
a1.text(tau / 2, -0.88, "$\\tau_m$", ha="center", va="top", fontsize=7, color=INK)
a1.text(1.5, -0.88, "frame $f$, $f{+}1$, $f{+}2$", ha="left", va="top", fontsize=6.6, color=GREY)
a1.text(0.0, 2.15, "true $\\nu > f_s/2$", ha="left", fontsize=6.6, color=RED, va="center")
a1.text(1.75, 2.15, "folded $\\beta$", ha="left", fontsize=6.6, color=GREY, va="center")
bg.text(xs[0] + ws[0] / 2, y0 + 0.30, "Same samples, two frequencies;\nphase correction at $\\beta$: $\\rho$ = 0.32",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.3)

# ---------- 2: method
cx = xs[1] + ws[1] / 2
bg.text(cx, y0 + 1.98, "$\\mathbf{Y}(\\beta)=\\sum_k\\, \\mathbf{D}(\\nu_k)\\,\\mathbf{J}\\,\\mathbf{X}_k$", ha="center", va="center", fontsize=9.5, color=INK)
bg.text(cx, y0 + 1.62, "$\\mathbf{D}(\\nu)=\\mathrm{diag}\\,(e^{\\,i2\\pi\\nu\\tau_m})$", ha="center", va="center", fontsize=8, color=INK)
bg.text(cx, y0 + 1.40, "true $\\nu$, not the folded $\\beta$", ha="center", va="center", fontsize=6.8, color=RED)
items = [("Single component", "phase derotation, cost of a\none-step reconstruction"),
         ("Colliding aliases", "joint solve, separable only\nwith serial timing"),
         ("Schedule $(\\delta, c, d)$", "estimated from the data\n(MAP, cross-validated)")]
for i, (h, b) in enumerate(items):
    yy = y0 + 1.08 - i * 0.40
    bg.add_patch(Rectangle((xs[1] + 0.12, yy - 0.15), 0.035, 0.30, fc=RED if i == 0 else ("#e08a1e" if i == 1 else BLUE), ec="none"))
    bg.text(xs[1] + 0.22, yy + 0.09, h, ha="left", va="center", fontsize=6.5, color=INK, fontweight="bold")
    bg.text(xs[1] + 0.22, yy - 0.06, b, ha="left", va="center", fontsize=5.9, color=GREY, linespacing=1.1)

# ---------- 3: data, with icons
def pig(cx, cy, s, col=PIG):
    bg.add_patch(Ellipse((cx, cy), 0.44 * s, 0.27 * s, fc=col, ec="none"))
    bg.add_patch(Circle((cx + 0.22 * s, cy + 0.04 * s), 0.12 * s, fc=col, ec="none"))
    bg.add_patch(Ellipse((cx + 0.33 * s, cy + 0.01 * s), 0.075 * s, 0.065 * s, fc=PIGL, ec="none"))
    bg.add_patch(Polygon([[cx + 0.15 * s, cy + 0.12 * s], [cx + 0.19 * s, cy + 0.22 * s], [cx + 0.25 * s, cy + 0.13 * s]], fc=PIGL, ec="none"))
    bg.add_patch(Circle((cx + 0.25 * s, cy + 0.07 * s), 0.013 * s, fc=INK, ec="none"))
    for dx in (-0.13, -0.05, 0.05, 0.13):
        bg.add_patch(FancyBboxPatch((cx + dx * s - 0.022 * s, cy - 0.19 * s), 0.044 * s, 0.09 * s,
                                    boxstyle="round,pad=0,rounding_size=0.01", fc=col, ec="none"))
    bg.add_patch(Arc((cx - 0.24 * s, cy + 0.04 * s), 0.07 * s, 0.07 * s, theta1=60, theta2=330, color=PIGL, lw=1.0))
def infant(cx, cy, s):
    bg.add_patch(Circle((cx, cy + 0.10 * s), 0.10 * s, fc=SKIN, ec="none"))
    bg.add_patch(FancyBboxPatch((cx - 0.10 * s, cy - 0.20 * s), 0.20 * s, 0.25 * s, boxstyle="round,pad=0,rounding_size=0.08", fc="#f3dd9a", ec="none"))
    bg.add_patch(Circle((cx - 0.035 * s, cy + 0.11 * s), 0.011 * s, fc=INK, ec="none")); bg.add_patch(Circle((cx + 0.035 * s, cy + 0.11 * s), 0.011 * s, fc=INK, ec="none"))
rows = [("pig", "HFOV pig", "44 fps, 10.3 Hz"),
        ("pig", "Pig, 2 conditions", "13 fps, 0.70 Hz series"),
        ("piglets", "2 piglets, 5 rec.", "13 fps, heart 2.6–3.7 Hz"),
        ("infant", "Infant", "13 fps, heart 2.5 Hz")]
for i, (kind, a, b) in enumerate(rows):
    yy = y0 + 2.05 - i * 0.40
    ix = xs[2] + 0.27
    if kind == "pig": pig(ix - 0.03, yy, 0.62)
    elif kind == "piglets": pig(ix - 0.10, yy, 0.42); pig(ix + 0.10, yy, 0.42, col=PIGL)
    else: infant(ix, yy, 0.62)
    bg.text(xs[2] + 0.56, yy + 0.07, a, ha="left", va="center", fontsize=6.5, color=INK, fontweight="bold")
    bg.text(xs[2] + 0.56, yy - 0.08, b, ha="left", va="center", fontsize=5.9, color=GREY)
bg.text(xs[2] + ws[2] / 2, y0 + 0.30, "9 records, 5 subjects;\nevery 2nd frame removed,\nfull-rate image = reference", ha="center", va="center", fontsize=6.2, color=INK, linespacing=1.2)

# ---------- 4: result: images (piglet 2 cardiac line) + paired plot over 17 components
rows4, cmap = FI.load("*_E7decim_E10_fig_g10-c005", "fig_decim_aligned")
panels = rows4[0]
show = [(0, "Reference"), (2, "Phase corr."), (5, "k-t")]
rho = {2: 0.83, 5: 1.00}
iw = 0.58
for j, (k, lab) in enumerate(show):
    ax = fig.add_axes([(xs[3] + 0.10 + j * 0.64) / W, (y0 + 1.47) / H, iw / W, iw / H])
    FI.draw_panel(ax, panels[k], cmap, "aligned")
    for ln in ax.lines: ln.set_linewidth(0.5)
    bg.text(xs[3] + 0.10 + j * 0.64 + iw / 2, y0 + 2.08, lab, ha="center", va="bottom", fontsize=6.4,
            color=RED if lab == "k-t" else INK, fontweight="bold")
    if k in rho:
        bg.text(xs[3] + 0.10 + j * 0.64 + iw / 2, y0 + 1.47, f"$\\rho$ = {rho[k]:.2f}", ha="center", va="top",
                fontsize=6.4, color=RED if lab == "k-t" else INK)
bg.text(xs[3] + 0.10 + iw / 2, y0 + 1.40, "piglet, cardiac\n3.5 Hz, 6.5 fps", ha="center", va="top", fontsize=5.6, color=GREY, linespacing=1.1)
# paired values from Table IV: (best baseline rho, k-t rho)
pairs = [(.73, .88), (.62, .89), (.47, .79), (.67, .75), (.50, .37), (.80, .95), (.82, .97), (.76, .93), (.65, .90), (.72, .55),
         (.78, .98), (.69, .88), (.72, .70), (.90, .99), (.96, 1.00), (.96, 1.00), (.97, 1.00)]
ap = fig.add_axes([(xs[3] + 0.33) / W, (y0 + 0.30) / H, 0.62 / W, 0.82 / H])
for b, k in pairs:
    up = k > b
    ap.plot([0, 1], [b, k], color=RED if up else GREY, lw=0.8 if up else 0.9, alpha=0.85 if up else 1, ls="-" if up else (0, (2, 1)))
    ap.plot([0, 1], [b, k], "o", ms=2.2, color=RED if up else GREY, mec="white", mew=0.3)
ap.set_xlim(-0.25, 1.25); ap.set_ylim(0.3, 1.02); ap.set_yticks([0.4, 0.7, 1.0])
ap.set_xticks([0, 1]); ap.set_xticklabels(["best base.", "k-t"], fontsize=5.8)
ap.tick_params(labelsize=6.0, length=2, width=0.5, pad=1.5); ap.set_ylabel("$\\rho$", fontsize=6.6, labelpad=1)
for s in ("top", "right"): ap.spines[s].set_visible(False)
for s in ("left", "bottom"): ap.spines[s].set_linewidth(0.6)
bg.text(xs[3] + 1.10, y0 + 1.00, "14 of 17\nimproved", ha="left", va="center", fontsize=7.4, color=RED, fontweight="bold", linespacing=1.15)
bg.text(xs[3] + 1.10, y0 + 0.66, "median $\\Delta\\rho$ +0.15\nWilcoxon $p$=0.003", ha="left", va="center", fontsize=5.8, color=INK, linespacing=1.2)
bg.text(xs[3] + 1.10, y0 + 0.30, "grey: 3 failures,\nweak lines near\nstrong slow lines", ha="left", va="center", fontsize=5.5, color=GREY, linespacing=1.15)

for ext in ("pdf", "svg"): fig.savefig(OUT / ("graphical_abstract." + ext))
fig.savefig(OUT / "graphical_abstract.png", dpi=360)
