"""Finalize T4: read the real OFF/ON rank-ablation JSONs and produce the trajectory
figure + summary. Usage: python fig4_gen.py off.json on.json"""
import json, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

off = json.loads(Path(sys.argv[1]).read_text())
on = json.loads(Path(sys.argv[2]).read_text())
OUT = Path("C:/Users/mobar/OneDrive/biorag/morpheus-paper/paper/figures")
OUT.mkdir(parents=True, exist_ok=True)

def series(rec):
    h = [r for r in rec["history"] if r["epoch"] >= 0]
    return [r["epoch"] for r in h], [r["biology_effective_rank"] for r in h]

ex, ey = series(off); nx, ny = series(on)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight", "axes.grid": True, "grid.alpha": 0.25})
fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.plot(ex, ey, "-o", ms=3, color="#c53030", label=f"OFF (decorr=0.0) — final {ey[-1]:.0f}")
ax.plot(nx, ny, "-o", ms=3, color="#2f855a", label=f"ON (decorr=0.04, F-R2) — final {ny[-1]:.0f}")
ax.set_xlabel("epoch"); ax.set_ylabel("biology-head effective rank (WSI, held-out)")
ax.set_title("T4 (real data): F-R2 decorrelation resists biology-head rank collapse")
ax.legend(frameon=False, loc="upper right")
fig.savefig(OUT / "fig4_ablation_trajectory.png"); plt.close(fig)

summary = {
    "off_final_rank": ey[-1], "on_final_rank": ny[-1], "delta": ny[-1] - ey[-1],
    "off_init": off["history"][0]["biology_effective_rank"], "on_init": on["history"][0]["biology_effective_rank"],
    "n_test": off.get("n_test"), "epochs": off.get("epochs"), "token_budget": off.get("token_budget"),
}
Path("C:/Users/mobar/OneDrive/biorag/morpheus-paper/paper/figures/t4_summary.json").write_text(json.dumps(summary, indent=2))
print("T4 SUMMARY:", json.dumps(summary))
print("OFF trajectory:", [round(v,1) for v in ey])
print("ON  trajectory:", [round(v,1) for v in ny])
print("wrote fig4_ablation_trajectory.png")
