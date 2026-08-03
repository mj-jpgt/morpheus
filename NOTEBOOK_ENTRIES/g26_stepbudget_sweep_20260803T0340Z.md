## 2026-08-03 03:40 UTC — CORRECTION: the variance floor delays the collapse rather than preventing it; lr 3e-3 is stable where lr 1e-3 is not

**Logged:** 2026-08-03 03:40 UTC. **How obtained:** `~/verify_varfix.py` on the Lambda A100 (`~/ws_d1`),
real cohort, one fixed 16-patient batch, hidden 512 / 4 layers / 8 heads, programme head 256, seed 42,
frozen memory queue, `variance` 0.01 beside `decorrelation` 0.04. Only steps and lr differ between arms.
Logs: `~/e0_run/collapse_fix.log`, `~/e0_run/g26_sweep.log`.

### Correction to the 02:10 UTC entry
That entry concluded "the collapse is gone" from an 800-step run. **That was wrong.** The collapse
still occurs at lr 1e-3; it simply arrives after step ~1,500, outside the window that was measured.
Treat the variance floor as necessary but not sufficient.

### Technical

| arm | steps | lr | trajectory | final InfoNCE | retrieval acc@1 | cross pos / neg cos | patient-to-patient cos |
|---|---|---|---|---:|---:|---|---:|
| D | 800 | 1e-3 | still descending at stop | 2.0875 | 0.125 | 0.9992 / 0.7291 | 0.7261 |
| E | 5000 | 1e-3 | 2.7731@500 → 2.4626@1000 → 2.7716@1500 → pinned from 2500 | 2.7726 | 0.062 | 0.9999 / 0.9999 | **1.0000** |
| F | 3000 | **3e-3** | 2.1184@500 → 2.1117@1000 → 2.1093@1500 → 2.1178@2000 → 2.1155@2500 | **2.1159** | **0.188** | 0.9983 / **0.6557** | **0.6480** |

Chance is ln(16) = 2.7726 for InfoNCE and 0.062 for retrieval. Reference: decorrelation deleted
entirely gives 2.0789 / 0.188 / 0.4946 at 800 steps, lr 1e-3.

**The learning-rate result is backwards from the usual stability story and is the strongest clue
available.** lr 1e-3 collapses by ~1,500 steps; lr 3e-3 does not collapse through 3,000 and its
patient-to-patient cosine *decreases* (0.7089 → 0.6480), i.e. the representation spreads out rather
than contracting. Arm F is a genuinely non-collapsing configuration, not a slower collapse. A
plausible mechanism — untested — is that the contrastive term escapes the collapse basin before the
decorrelation attractor captures it, and that at the lower rate decorrelation wins the race.

**Arm F is still ~20× above the `<= 0.10` gate criterion.** It is a starting point, not a solution.

Attainability arithmetic, so the gate is not assumed impossible: at temperature 0.07, positives at
cosine 0.999 against negatives near 0.5 give a loss of ≈0.06. The blocker is that negatives sit at
0.66, not that the criterion is unreachable.

### In plain terms
Adding the anti-collapse rule back bought time rather than a cure — at the slow learning rate the
model still ends up making every patient look identical, just later than we first measured. Training
*faster* avoids it entirely, which is the opposite of what one would expect and is the best hint we
have about what is really going on.

### Meaning for the claim
- D1 stays blocked; three defects now, two fixed, one open.
- The G2.6 threshold has not been touched and, by the arithmetic above, does not need to be.
- Any D1 launch configuration must be shown stable over the FULL step budget, not a window — this
  correction exists because a window was mistaken for a fix.

### Files / commits
`v2/losses.py`, `v2/training.py`, `v2/tests/test_programme_free.py`; `~/verify_varfix.py`.
Supersedes the conclusion of `g26_variance_floor_fix_20260803T0210Z.md`.
