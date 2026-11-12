#!/usr/bin/env python3
"""Figure 4: the kidney-function pathway from AITD liability to erythroid traits."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
med = pd.read_csv(os.path.join(BASE, "results", "mediation.tsv"), sep="\t")
k = med[med.mediator == "k"].set_index(["exposure", "outcome"])
LN2 = np.log(2)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

C_AIH, C_GD = "#B03A2E", "#2E5E8C"
C_TBIL, C_MRV = "#C8912B", "#4A7C59"

fig = plt.figure(figsize=(7.25, 3.15))
ax1 = fig.add_axes([0.005, 0.02, 0.60, 0.92])
ax2 = fig.add_axes([0.695, 0.17, 0.30, 0.71])

# ---------------------------------------------------------------- panel A
ax1.set_xlim(0, 10); ax1.set_ylim(0, 6); ax1.axis("off")


def box(x, y, w, h, text, fc="#FFFFFF", ec="#333333", fs=6.4, bold=False):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                       linewidth=0.8, edgecolor=ec, facecolor=fc, zorder=2)
    ax1.add_patch(b)
    ax1.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
             fontweight="bold" if bold else "normal", zorder=3, linespacing=1.4)


def arrow(x0, y0, x1, y1, color="#333333", lw=0.9, rad=0.0):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
                        linewidth=lw, color=color, zorder=1,
                        connectionstyle=f"arc3,rad={rad}", shrinkA=1, shrinkB=1)
    ax1.add_patch(a)


b1 = k.loc[("AIH", "TBIL"), "b1"] * LN2
se1 = k.loc[("AIH", "TBIL"), "se1"] * LN2

box(0.10, 2.05, 2.75, 1.70,
    "Genetic liability to\nautoimmune\nhypothyroidism\n(per doubling of odds)",
    fc="#EAF0F6")
box(4.25, 2.25, 2.15, 1.30, "Serum creatinine\n(per SD)\nhigher = worse", fc="#F5EDE6")

arrow(2.88, 2.90, 4.22, 2.90, lw=1.1)
ax1.text(3.55, 3.80,
         f"β₁ = +{b1:.3f}\n({b1 - 1.96 * se1:.3f} to {b1 + 1.96 * se1:.3f})",
         ha="center", va="bottom", fontsize=5.9, linespacing=1.3)

rows = [("TBIL", "Total bilirubin", 4.15, C_TBIL), ("MRV", "Mean reticulocyte\nvolume", 0.45, C_MRV)]
for oc, label, yy, col in rows:
    d = k.loc[("AIH", oc)]
    prop = d["prop"] * 100
    box(7.55, yy, 2.30, 1.55,
        f"{label}\nβ₂ = {d['b2']:+.3f}\n{prop:.0f}% mediated",
        fc="#FFFFFF", ec=col, fs=6.2)
    arrow(6.42, 2.90, 7.53, yy + 0.78, color=col, lw=1.5, rad=0.15 if yy < 3 else -0.15)

ax1.text(0.10, 5.78, "A", fontsize=9, fontweight="bold", va="top")
ax1.text(0.52, 5.44, "Kidney-function pathway (autoimmune hypothyroidism)",
         fontsize=7.2, fontweight="bold", va="top")
def _ind(oc):
    d = k.loc[("AIH", oc)]
    return (d["ind"] * LN2, d["ind_lo"] * LN2, d["ind_hi"] * LN2)


_t, _tl, _th = _ind("TBIL")
_m, _ml, _mh = _ind("MRV")
ax1.text(0.10, 0.05,
         f"Indirect effect: total bilirubin {_t:+.4f} SD (95% CI {_tl:.4f} to {_th:.4f});\n"
         f"mean reticulocyte volume {_m:+.4f} SD ({_ml:.4f} to {_mh:.4f}).",
         fontsize=6.0, va="bottom", color="#333333", linespacing=1.4)

# ---------------------------------------------------------------- panel B
ORDER = ["TBIL", "MRV", "MCV", "DBIL", "RDW"]
LAB = ["Total\nbilirubin", "Mean retic.\nvolume", "MCV", "Direct\nbilirubin", "RDW"]
x = np.arange(len(ORDER))

for off, ex, col, name in [(-0.19, "AIH", C_AIH, "Autoimmune hypothyroidism"),
                           (0.19, "GD", C_GD, "Graves' disease")]:
    v = np.array([k.loc[(ex, o), "prop"] * 100 for o in ORDER])
    lo = np.array([k.loc[(ex, o), "prop"] - k.loc[(ex, o), "prop_lo"] for o in ORDER]) * 100
    hi = np.array([k.loc[(ex, o), "prop_hi"] - k.loc[(ex, o), "prop"] for o in ORDER]) * 100
    ax2.bar(x + off, v, width=0.34, color=col, edgecolor="none", label=name, zorder=2)
    ax2.errorbar(x + off, v, yerr=[lo, hi], fmt="none", ecolor="#333333",
                 elinewidth=0.7, capsize=1.8, zorder=3)

ax2.axhline(0, color="#333333", lw=0.7)
ax2.set_xticks(x); ax2.set_xticklabels(LAB, fontsize=6.4, linespacing=1.2)
ax2.set_ylabel("% of total effect mediated\nby kidney function", fontsize=6.8)
ax2.set_ylim(-62, 90)
ax2.tick_params(axis="y", labelsize=6.2)
ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
ax2.legend(fontsize=6.0, frameon=False, loc="upper right", handlelength=1.0,
           borderaxespad=0.1)
ax2.set_title("B  Mediated proportion by outcome", fontsize=7.4, fontweight="bold", loc="left")

for _ext in ("png", "pdf"):
    plt.savefig(os.path.join(BASE, "figures", f"fig4_pathway.{_ext}"), dpi=400,
                facecolor="white", bbox_inches="tight", pad_inches=0.02)
print("saved figures/fig4_pathway.png and .pdf")
