"""Generate MORPHEUS rank-collapse paper figures (T1-T3) from the verified numbers."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUT = Path("C:/Users/mobar/OneDrive/biorag/morpheus-paper/paper/figures")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight", "axes.grid": True,
                     "grid.alpha": 0.25, "grid.linestyle": "-"})
BLUE, RED, GREY, GREEN = "#2b6cb0", "#c53030", "#718096", "#2f855a"

# ---- Fig 1: dual-head effective-rank spectrum (T1) ----
states = ["wsi\nidentity", "rna\nidentity", "full\nidentity", "wsi\nbiology", "rna\nbiology", "full\nbiology", "full\npatient"]
ranks = [84.3, 37.5, 55.8, 6.0, 4.4, 5.3, 8.0]
colors = [BLUE, BLUE, BLUE, RED, RED, RED, GREY]
fig, ax = plt.subplots(figsize=(7.2, 3.8))
bars = ax.bar(range(len(states)), ranks, color=colors)
for i, v in enumerate(ranks):
    ax.text(i, v + 2, f"{v:.1f}", ha="center", fontsize=9)
ax.set_xticks(range(len(states))); ax.set_xticklabels(states)
ax.set_ylabel("effective rank (Roy–Vetterli)  /  256 dims"); ax.set_ylim(0, 95)
ax.set_title("Dual-head rank fingerprint: biology heads collapse (~5), identity stays healthy (~84)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE, label="identity heads"), Patch(color=RED, label="biology heads"),
                   Patch(color=GREY, label="fused patient")], loc="upper right", frameon=False, fontsize=9)
fig.savefig(OUT / "fig1_rank_spectrum.png"); plt.close(fig)

# ---- Fig 2: confound decomposition ladder (T2) ----
methods = ["MLP-CLIP", "SigLIP", "MORPHEUS-v2\n(anchored)", "MORPHEUS-v2\n(no-anchor)"]
global_r = [0.348, 0.349, 0.327, 0.338]
within_r = [0.188, 0.193, 0.166, 0.185]
random_within = 0.154  # random gene-set control (within-cancer, MLP-CLIP)
x = np.arange(len(methods)); w = 0.38
fig, ax = plt.subplots(figsize=(7.2, 4.0))
ax.bar(x - w/2, global_r, w, color=GREY, label="pooled (cross-cancer)")
ax.bar(x + w/2, within_r, w, color=BLUE, label="within-cancer")
ax.axhline(random_within, ls="--", c=RED, lw=1.3, label="random gene-set null (within-cancer)")
for i in range(len(methods)):
    drop = 100 * (global_r[i] - within_r[i]) / global_r[i]
    ax.annotate(f"-{drop:.0f}%", (x[i], within_r[i]), (x[i], global_r[i] + 0.012),
                ha="center", fontsize=8.5, color=RED,
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
ax.set_xticks(x); ax.set_xticklabels(methods); ax.set_ylabel("Pearson r (Hallmark prompting)")
ax.set_ylim(0, 0.40); ax.legend(loc="upper right", frameon=False, fontsize=9)
ax.set_title("~46–49% of pooled Pearson is cross-cancer cohort structure (all methods)")
fig.savefig(OUT / "fig2_confound_ladder.png"); plt.close(fig)

# ---- Fig 3: method-invariant control-adjusted specificity (T3) ----
spec_within = [0.068, 0.067, 0.068, 0.069]  # real - random, within-cancer
allm = methods + ["MLP-CLIP\nhard-neg", "MORPHEUS-v1"]
spec_all = spec_within + [0.066, 0.070]
fig, ax = plt.subplots(figsize=(7.2, 3.8))
ax.bar(range(len(allm)), spec_all, color=GREEN)
mean = np.mean(spec_all)
ax.axhline(mean, ls="--", c="k", lw=1.0)
ax.text(len(allm) - 0.5, mean + 0.002, f"mean ≈ +{mean:.3f}", ha="right", fontsize=9)
for i, v in enumerate(spec_all):
    ax.text(i, v + 0.001, f"+{v:.3f}", ha="center", fontsize=8.5)
ax.set_xticks(range(len(allm))); ax.set_xticklabels(allm, fontsize=8.5)
ax.set_ylabel("control-adjusted within-cancer\nspecificity (Δ Pearson)")
ax.set_ylim(0, 0.09)
ax.set_title("Genuine biology signal is ~+0.07 and method-invariant (incl. baseline)")
fig.savefig(OUT / "fig3_method_invariance.png"); plt.close(fig)

print("wrote:", *[p.name for p in sorted(OUT.glob('*.png'))])
