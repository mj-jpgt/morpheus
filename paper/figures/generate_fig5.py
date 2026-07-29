"""Finalize the multi-seed T4 sweep: read s{42,43,44}_{off,on}.json, summarize
rank recovery + specificity across seeds, plot fig5, print paper-update text.
Usage: python finalize_multiseed.py <dir-with-jsons>"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = Path(sys.argv[1])
OUT = Path("C:/Users/mobar/OneDrive/biorag/morpheus-paper/paper/figures")
SEEDS = [42, 43, 44]

def load(seed, arm):
    p = D / f"s{seed}_{arm}.json"
    if not p.exists(): return None
    return json.loads(p.read_text())

rows = []
for s in SEEDS:
    off, on = load(s, "off"), load(s, "on")
    if off is None or on is None:
        print(f"MISSING seed {s}: off={off is not None} on={on is not None}"); continue
    r = {
        "seed": s,
        "off_rank": off["final_biology_effective_rank"],
        "on_rank": on["final_biology_effective_rank"],
        "off_spec": (off.get("prompting") or {}).get("prompting_within_cancer_pearson"),
        "on_spec": (on.get("prompting") or {}).get("prompting_within_cancer_pearson"),
        "off_pool": (off.get("prompting") or {}).get("prompting_pooled_pearson"),
        "on_pool": (on.get("prompting") or {}).get("prompting_pooled_pearson"),
    }
    r["d_rank"] = r["on_rank"] - r["off_rank"]
    rows.append(r)

def ms(key):
    v = np.array([r[key] for r in rows if r[key] is not None], float)
    return (np.mean(v), np.std(v)) if len(v) else (float("nan"), float("nan"))

print("\n=== PER-SEED ===")
for r in rows:
    print(f" seed {r['seed']}: rank off={r['off_rank']:.1f} on={r['on_rank']:.1f} (Δ+{r['d_rank']:.1f}) | "
          f"within-cancer spec off={r['off_spec']} on={r['on_spec']}")
mo, so_ = ms("off_rank"); mn, sn = ms("on_rank"); md, sd = ms("d_rank")
print(f"\n=== ACROSS {len(rows)} SEEDS ===")
print(f" biology rank  OFF {mo:.1f}±{so_:.1f}   ON {mn:.1f}±{sn:.1f}   Δ +{md:.1f}±{sd:.1f}")
osm, oss = ms("off_spec"); nsm, nss = ms("on_spec")
print(f" within-cancer specificity  OFF {osm:.3f}±{oss:.3f}   ON {nsm:.3f}±{nss:.3f}")

# fig5: per-seed OFF vs ON rank
fig, ax = plt.subplots(figsize=(6.4, 4.0))
x = np.arange(len(rows)); w = 0.38
ax.bar(x - w/2, [r["off_rank"] for r in rows], w, color="#c53030", label="baseline (variance floor)")
ax.bar(x + w/2, [r["on_rank"] for r in rows], w, color="#2f855a", label="F-R2 (decorrelation)")
for i, r in enumerate(rows):
    ax.text(i - w/2, r["off_rank"]+1, f"{r['off_rank']:.0f}", ha="center", fontsize=8)
    ax.text(i + w/2, r["on_rank"]+1, f"{r['on_rank']:.0f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([f"seed {r['seed']}" for r in rows])
ax.set_ylabel("biology-head effective rank (WSI, held-out)")
ax.set_title(f"Multi-seed: F-R2 recovers biology rank (+{md:.0f}±{sd:.0f}) across seeds")
ax.legend(frameon=False, loc="upper left", fontsize=9)
ax.spines[["top","right"]].set_visible(False); ax.grid(axis="y", alpha=0.25)
fig.savefig(OUT / "fig5_multiseed_rank.png", dpi=150, bbox_inches="tight"); plt.close(fig)
Path(OUT / "t4_multiseed_summary.json").write_text(json.dumps({
    "seeds": rows, "off_rank_mean": mo, "off_rank_std": so_, "on_rank_mean": mn, "on_rank_std": sn,
    "d_rank_mean": md, "d_rank_std": sd, "off_spec_mean": osm, "on_spec_mean": nsm}, indent=2))
print("\nwrote fig5_multiseed_rank.png + t4_multiseed_summary.json")
