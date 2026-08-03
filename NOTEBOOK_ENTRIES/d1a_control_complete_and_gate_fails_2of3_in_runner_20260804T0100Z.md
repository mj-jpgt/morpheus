## 2026-08-04 01:00 UTC — D1-A closes as the control: collapse confirmed at epoch 40, and the gate passes 3/3 in the harness but 1/3 in the runner

**Logged:** 2026-08-04 01:00 UTC. **How obtained:** `v2/research/rebase/d1_geometry_probe.py` on the
A100 (`150.136.45.194`), final epoch-40 checkpoints, 282 held-out test patients;
`~/e0_run/d1_v1_launch.log` for the gate outcomes.

### Technical

**A5 — effective rank at epoch 40.** Centred participation ratio of singular values of `z_biology`
(WSI view) on held-out test patients:

| arm | epoch | centred eff-rank | hard rank @1e-3 | RNA–RNA cos | feat-std |
|---|---|---:|---:|---:|---:|
| `programme_only` seed 42 | 39 | **9.81** | 136 | 0.258 | 0.0056 |
| `programme_only` seed 43 | 39 | **10.47** | 158 | 0.265 | 0.0063 |
| `programme_free` seed 42 | 39 | **1.71** | 11 | 0.986 | 0.0156 |

The supervised arm **improves with training** (7.38 → 9.81 between epochs 21 and 40) while the
contrastive arm **stays collapsed** (1.76 → 1.71). The final gap is **5.7–6.1×**, and the contrastive
arm's RNA-view states sit at mutual cosine 0.986 — very nearly one vector across 282 patients.

This is what D1-A was kept running for. The defect is now documented at full training duration rather
than inferred from a mid-run projection, and the two `programme_only` seeds agree (9.81 / 10.47),
so the comparison arm is stable.

**D1-A's final disposition:**

| run | outcome |
|---|---|
| `d1_p_seed42` | 40 epochs ✓ |
| `d1_f_seed42` | 40 epochs ✓ — **the control pair** |
| `d1_p_seed43` | 40 epochs ✓ |
| `d1_f_seed43` | G2.6 failed, contrastive **0.50883** |
| `d1_p_seed44` | 27 epochs when the pipeline unwound |
| `d1_f_seed44` | G2.6 failed, contrastive **2.14122** |

`run_d1` raises on the first non-zero return code, so no exports, CALIBRA or bootstrap were produced.
Seed 42's complete pair is intact and is sufficient for the control's purpose.

**The harness-versus-runner gap is far larger than one data point suggested.** The same gate function,
same seeds, same 2,400-step budget:

| seed | standalone harness | inside the runner |
|---|---:|---:|
| 42 | 0.01871 | **passed** |
| 43 | 0.01206 | **0.50883** ✗ |
| 44 | 0.05666 | **2.14122** ✗ |

**3/3 pass in the harness; 1/3 in the runner.** Seed 44's in-runner value of 2.14 is not marginal — it
is close to the chance value ln(16) = 2.7726. The harness did not merely give optimistic numbers; it
inverted the verdict on two of three seeds, and I used it to decide to launch. This strengthens
instance 4 of `paper/LIVENESS_GATE_DESIGN.md` from an anecdote to a rate.

It also, incidentally, means the repaired gate is **doing its job**: it refused to train two arms whose
objective we now independently know collapses to effective rank ~1.7. The gate was right and the
harness was wrong.

### In plain terms

The control run finished what it was kept alive to do. At the end of full training, the supervised arm
describes a held-out patient with about ten independent numbers and is still improving; the
unsupervised arm uses fewer than two and has not moved since epoch 21. That is the defect, measured at
full duration rather than guessed at from halfway.

Separately, the pre-flight check I used to decide the run was safe to launch gave a clean pass on all
three repeats, and the real check inside the training job then rejected two of them — one of those
almost at the "random guessing" score. The stand-in did not just flatter the result; it reversed it.

### Meaning for the claim

The D1 collapse is now a documented result at full training duration, with a stable two-seed
comparison arm, rather than a mid-run projection. Nothing about programme supervision may be concluded
from it — the contrastive arm never trained, so the comparison measures a defect, not an ablation.

D1-B is still correctly held pending the momentum durability runs.

### Files / commits

- `v2/research/rebase/d1_geometry_probe.py`
- `~/e0_run/d1_v1/d1_{p,f}_seed42/last.pt` (epoch 40), `d1_p_seed43/last.pt`
- `~/e0_run/d1_v1_launch.log`
- Strengthens: `paper/LIVENESS_GATE_DESIGN.md` instance 4
