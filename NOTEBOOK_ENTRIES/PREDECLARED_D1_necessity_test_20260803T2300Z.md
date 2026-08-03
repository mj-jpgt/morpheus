## 2026-08-03 23:00 UTC — PREDECLARATION: how D1's rank-vs-channel comparison will be read, all four outcomes, before any channel number is seen

**Logged:** 2026-08-03 23:00 UTC. **Written and committed before opening `d1_audit.log`, before any
CALIBRA channel value has been read, and before the paired bootstrap has run.** Effective ranks below
are already measured and are stated here so the prediction is anchored to real numbers; **the channel
side is entirely unseen.**

**This predeclaration exists because the author has told me in advance which result supports the
paper's claim.** That is the correct thing for them to have done and it is exactly the circumstance in
which a prior commitment to the reading is worth more than the reading itself.

### What is being tested

P2's broad claim is that *effective rank does not track information content*. RankMe's defence is that
high rank is **"a necessary (but not sufficient) condition"**. Under that hedge, every *high rank +
low information* observation is already predicted and is not a counterexample — and LiDAR
(arXiv:2312.04000) has in any case already published that direction. The result that breaks necessity
is the reverse: **low rank carrying high information.**

D1 is the test: `programme_free` has lower effective rank than `programme_only`; does it carry a
comparable or better molecular channel?

### Correction to the premise, made before the test

Two things must be fixed or the test is misreported.

**1. The rank gap is ~2×, not ~9×.** The "12 vs 111" figures are the **in-run tripwire**: participation
ratio, on *training* batches, at *step 200*. The quantity the channel is computed on is different —
*held-out* test patients, at *epoch 40*, on the exported artifact. Measured there just now:

| seed | `programme_free` | `programme_only` | ratio |
|---|---:|---:|---:|
| 42 | 11.67 | 22.18 | 1.9× |
| 43 | 6.60 | 18.85 | 2.9× |
| 44 | 5.86 | 9.34 | 1.6× |

*(canonical definition, 282 held-out patients, epoch 39)*

So `programme_free` is **1.6–2.9× lower**, not 9×. The necessity test is still meaningful — the arms
are separated and consistently ordered — but "9× lower rank" would be a false statement in the paper.

**2. My earlier numbers use a different definition from the canonical one.** Two implementations exist
and they are not the same function:

| definition | formula | used by |
|---|---|---|
| **canonical** — `v2/calibra/spectral.py`, Roy–Vetterli, **and RankMe's** | `exp(−Σ p log p)`, `p = σ/Σσ` | CALIBRA, and this predeclaration |
| participation ratio | `(Σσ)² / Σσ²` | `d1_geometry_probe.py`, `d1_audit.py`, the in-run tripwire |

Measured on the same states they differ by 1.5–2.2×, and the ratio is **not constant** (0.46–0.69), so
they are not interchangeable even up to a scale factor. **All rank numbers I report for the paper will
use the canonical definition**, because a paper about RankMe must quote RankMe's. The tripwire keeps
the participation ratio, which is fine — it is an operational abort criterion, not a paper number —
but that is now stated rather than assumed. The three-implementation reconciliation is another agent's;
this entry records which one each of my numbers uses so they can reconcile without re-deriving it.

### The four outcomes, and how each will be read

Let Δ = channel(`programme_free`) − channel(`programme_only`), per seed, with its bootstrap CI.

| # | outcome | reading |
|---|---|---|
| **O1** | Δ ≥ 0, or CI includes 0, on **2 of 3 or 3 of 3** seeds | **Necessity violated.** Lower rank, equal-or-better channel. The broad claim survives its strongest objection and this is P2's headline. |
| **O2** | Δ decisively negative (CI excludes 0) on **3 of 3** seeds | **Rank vindicated on this pair.** Goes in the paper as a limitation, with the same prominence a confirmation would have had. The broad claim must then be narrowed or defended on other evidence. |
| **O3** | Δ mixed — decisively negative on some seeds, non-negative on others | **Rank is unreliable as a selector**, which is weaker than O1 but still incompatible with a *necessary* condition used for model selection. Report the per-seed split, not a pooled mean. |
| **O4** | ranks not actually separated once channel-matched, or CIs too wide to distinguish anything | **The pair does not test necessity.** Report as uninformative. Do not convert a null into support. |

**The operational headline, predeclared:** *how often would a rank-based selection rule have chosen
the arm with the larger channel?* The rule picks the higher-rank arm, which is `programme_only` on all
three seeds. So the fraction is the number of seeds where `programme_only` also has the larger channel,
out of 3.

- **3/3** → rank is a reliable selector here; supports RankMe on this pair.
- **0/3** → rank is anti-correlated with the channel; strongest form of the claim.
- **1/3 or 2/3** → rank is no better than a coin flip on this pair. With n=3 this is weak evidence and
  will be reported as weak, not as a headline.

**Stated in advance so it cannot be chosen afterwards:** with three seeds, 1/3 and 2/3 are not
publishable as a rate. Only 0/3 and 3/3 are worth a sentence, and even those are three paired
comparisons, not a benchmark.

### What would make me distrust a favourable result

- If `programme_free`'s channel is high **only** on seed 42, whose rank (11.67) is nearly double the
  other two arms' (6.60, 5.86), the effect tracks rank *within* the arm and argues for necessity, not
  against it.
- If the channel difference is within the noise CALIBRA's own controls establish, it is not a
  difference.
- If `random_control` shows any channel, none of this is interpretable, per the standing audit rule.

### Files / commits

- Ranks above: `/tmp/rank_defs.py` on the box, D1-B checkpoints at epoch 39, 282 held-out patients
- `v2/calibra/spectral.py` (canonical), `v2/research/rebase/d1_geometry_probe.py` (participation ratio)
- Prior art: `NOTEBOOK_ENTRIES/p2_rank_draft_20260803T2134Z.md` — RankMe's hedge, and LiDAR having
  already published the high-rank/low-information direction
- Same discipline previously applied: `PREDECLARED_turnover_criterion` (falsified),
  `PREDECLARED_learning_rate_test` (pending)
