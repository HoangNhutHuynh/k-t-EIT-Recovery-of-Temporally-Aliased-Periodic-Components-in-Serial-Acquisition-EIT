"""Graphical abstract of the k-t EIT manuscript for IEEE Xplore (file name gagraphic).
IEEE spec: 660 x 295 px, sans-serif text >= 16 px, labels only, contrast >= 4.5:1.
Layout is written in pixel units of the 660 x 295 canvas; 1 px = 0.8 pt at 90 dpi.
The numbers in panel 3 are those of the frame-removal table of the manuscript.

Part of k-t EIT. Author: Hoang Nhut Huynh. License: MIT.
Run from the repository root:  python figures/graphical_abstract.py
Uses the image data of the newest '*_E7decim_E10_fig_g10-c005' run (see run_export_fig_data.m);
writes gagraphic.png (660 x 295), gagraphic_300dpi.png and gagraphic_review.pdf next to this script.
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle, FancyBboxPatch, Polygon, Rectangle, FancyArrowPatch, Arc
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_images as FI

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "mathtext.fontset": "dejavusans"})
PX = 0.8                                   # pt per px at 90 dpi
def fs(px): return px * PX
INK, GREY, RED, BLUE, EDGE = "#1f2937", "#4b5563", "#b42318", "#1d4ed8", "#c9cfd8"
PIG, PIGL, SKIN = "#f4b6a6", "#e79c88", "#f5d1ad"
Wpx, Hpx = 660, 295
fig = plt.figure(figsize=(Wpx / 90, Hpx / 90))
bg = fig.add_axes([0, 0, 1, 1]); bg.set_xlim(0, Wpx); bg.set_ylim(0, Hpx); bg.axis("off")
def ax_px(x, y, w, h): return fig.add_axes([x / Wpx, y / Hpx, w / Wpx, h / Hpx])

# ---- three panels
xs, ws = [6, 226, 442], [206, 202, 212]
y0, hh = 6, 283
heads = ["Serial EIT: aliasing", "k-t forward model", "Result"]
for k, (x, w, t) in enumerate(zip(xs, ws, heads)):
    bg.add_patch(FancyBboxPatch((x, y0), w, hh, boxstyle="round,pad=0,rounding_size=8", fc="white", ec=EDGE, lw=1.0))
    bg.text(x + w / 2, y0 + hh - 16, t, ha="center", va="center", fontsize=fs(18), fontweight="bold", color=INK)
    if k < 2:
        bg.add_patch(FancyArrowPatch((x + w - 2, y0 + hh / 2), (xs[k + 1] + 2, y0 + hh / 2),
                                     arrowstyle="-|>", mutation_scale=14, color="#6b7280", lw=1.5, zorder=5))

# ---- 1: one measurement sampled once per frame; true vs folded frequency
a1 = ax_px(14, 70, 190, 168); a1.axis("off"); a1.set_xlim(-0.05, 3.05); a1.set_ylim(-1.0, 2.6)
t = np.linspace(0, 3, 600); tau = 0.40
a1.plot(t, 1.15 + 0.62 * np.cos(2 * np.pi * 1.25 * t), color=RED, lw=1.8)
a1.plot(t, 1.15 + 0.62 * np.cos(2 * np.pi * 0.25 * t + 2 * np.pi * tau), color=GREY, lw=1.6, ls=(0, (3, 1.6)))
ts = np.arange(3) + tau
a1.plot(ts, 1.15 + 0.62 * np.cos(2 * np.pi * 1.25 * ts), "o", color=INK, ms=6, zorder=5)
for f in range(3):
    a1.add_patch(Rectangle((f, -0.55), 1, 0.55, fc="#eef1f5", ec="#9aa3b2", lw=0.8))
    for s in range(8):
        a1.add_patch(Rectangle((f + 0.03 + s * 0.115, -0.49), 0.085, 0.43, fc="#1f4f8f" if s == 3 else "#9db9e0", ec="none"))
    a1.plot([f + tau, f + tau], [0.0, 1.15 + 0.62 * np.cos(2 * np.pi * 1.25 * (f + tau))], color="#9aa3b2", lw=0.8, ls=":")
a1.text(0.0, 2.35, "true $\\nu$", ha="left", va="center", fontsize=fs(16), color=RED, fontweight="bold")
a1.text(3.0, 2.35, "folded $\\beta$", ha="right", va="center", fontsize=fs(16), color=GREY, fontweight="bold")
bg.text(xs[0] + ws[0] / 2, 52, "serial times $\\tau_m$", ha="center", va="center", fontsize=fs(16), color=INK)
bg.text(xs[0] + ws[0] / 2, 22, "> $f_s/2$: wrong image", ha="center", va="center", fontsize=fs(16), color=INK, fontweight="bold")

# ---- 2: model + data
cx = xs[1] + ws[1] / 2
bg.text(cx, 222, "$\\mathbf{D}(\\nu)=\\mathrm{diag}(e^{\\,i2\\pi\\nu\\tau_m})$", ha="center", va="center", fontsize=fs(19), color=INK)
bg.text(cx, 190, "phase at the true $\\nu$", ha="center", va="center", fontsize=fs(16), color=RED, fontweight="bold")
bg.text(cx, 162, "joint solve of collisions", ha="center", va="center", fontsize=fs(16), color=INK)
bg.plot([xs[1] + 18, xs[1] + ws[1] - 18], [140, 140], color=EDGE, lw=1.0)
def pig(cx, cy, s, col=PIG):
    bg.add_patch(Ellipse((cx, cy), 44 * s, 27 * s, fc=col, ec="none"))
    bg.add_patch(Circle((cx + 22 * s, cy + 4 * s), 12 * s, fc=col, ec="none"))
    bg.add_patch(Ellipse((cx + 33 * s, cy + 1 * s), 7.5 * s, 6.5 * s, fc=PIGL, ec="none"))
    bg.add_patch(Polygon([[cx + 15 * s, cy + 12 * s], [cx + 19 * s, cy + 22 * s], [cx + 25 * s, cy + 13 * s]], fc=PIGL, ec="none"))
    bg.add_patch(Circle((cx + 25 * s, cy + 7 * s), 1.4 * s, fc=INK, ec="none"))
    for dx in (-13, -5, 5, 13):
        bg.add_patch(Rectangle((cx + (dx - 2.2) * s, cy - 19 * s), 4.4 * s, 9 * s, fc=col, ec="none"))
    bg.add_patch(Arc((cx - 24 * s, cy + 4 * s), 7 * s, 7 * s, theta1=60, theta2=330, color=PIGL, lw=1.2))
def infant(cx, cy, s):
    bg.add_patch(Circle((cx, cy + 10 * s), 10 * s, fc=SKIN, ec="none"))
    bg.add_patch(FancyBboxPatch((cx - 10 * s, cy - 20 * s), 20 * s, 25 * s, boxstyle="round,pad=0,rounding_size=7", fc="#f3dd9a", ec="none"))
    for d in (-3.5, 3.5): bg.add_patch(Circle((cx + d * s, cy + 11 * s), 1.2 * s, fc=INK, ec="none"))
pig(xs[1] + 36, 100, 0.75)
pig(xs[1] + 88, 95, 0.5); pig(xs[1] + 122, 95, 0.5, col=PIGL); pig(xs[1] + 156, 95, 0.5)
infant(xs[1] + 188, 100, 0.85)
bg.text(cx, 58, "HFOV pig, 3 piglets,\ninfant", ha="center", va="center", fontsize=fs(16), color=INK, linespacing=1.15)
bg.text(cx, 20, "9 records, 5 subjects", ha="center", va="center", fontsize=fs(16), color=INK, fontweight="bold")

# ---- 3: images of an aliased cardiac line + summary
rows4, cmap = FI.load("*_E7decim_E10_fig_g10-c005", "fig_decim_aligned")
panels = rows4[0]
show = [(0, "Ref.", None), (2, "Phase c.", 0.83), (5, "k-t", 1.00)]
iw = 64
for j, (k, lab, r) in enumerate(show):
    x = xs[2] + 8 + j * 68
    ax = ax_px(x, 148, iw, iw)
    FI.draw_panel(ax, panels[k], cmap, "aligned")
    for ln in ax.lines: ln.set_linewidth(0.8)
    col = RED if lab == "k-t" else INK
    bg.text(x + iw / 2, 226, lab, ha="center", va="center", fontsize=fs(16), color=col, fontweight="bold")
    if r is not None:
        bg.text(x + iw / 2, 134, f"{r:.2f}", ha="center", va="center", fontsize=fs(16), color=col, fontweight="bold")
bg.text(xs[2] + 8 + iw / 2, 134, "$\\rho$:", ha="center", va="center", fontsize=fs(16), color=INK)
bg.plot([xs[2] + 18, xs[2] + ws[2] - 18], [114, 114], color=EDGE, lw=1.0)
bg.text(xs[2] + ws[2] / 2, 88, "best in 14 of 17", ha="center", va="center", fontsize=fs(19), color=RED, fontweight="bold")
bg.text(xs[2] + ws[2] / 2, 56, "aliased components", ha="center", va="center", fontsize=fs(16), color=INK)
bg.text(xs[2] + ws[2] / 2, 22, "Wilcoxon p = 0.003", ha="center", va="center", fontsize=fs(16), color=INK)

if __name__ == "__main__":
    from pathlib import Path
    OUT = Path(__file__).resolve().parent
    fig.savefig(OUT / "gagraphic_review.pdf")
    fig.savefig(OUT / "gagraphic_300dpi.png", dpi=300)
    fig.savefig(OUT / "gagraphic.png", dpi=90)      # 660 x 295 px
