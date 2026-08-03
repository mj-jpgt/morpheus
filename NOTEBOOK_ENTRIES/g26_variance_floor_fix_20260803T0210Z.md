## 2026-08-03 02:10 UTC — The variance floor eliminates the programme_free collapse; G2.6 still not passing

**Logged:** 2026-08-03 02:10 UTC. **How obtained:** `~/verify_varfix.py` on the Lambda A100 (`~/ws_d1`),
real cohort, one fixed 16-patient batch, hidden 512 / 4 layers / 8 heads, programme head 256, seed 42,
800 steps, lr 1e-3, frozen memory queue. Identical to the arm that collapsed except that `variance`
is now active beside `decorrelation`.

### Technical

| metric | full schedule (before) | **+ variance floor** | decorrelation deleted (reference) |
|---|---|---|---|
| in-batch InfoNCE (chance 2.7726) | 2.7734 | **2.0875** | 2.0789 |
| retrieval acc@1 (chance 0.062) | 0.000 | **0.125** | 0.188 |
| cross-modal positive cos | 0.9959 | 0.9992 | 0.9988 |
| cross-modal negative cos | 0.9960 | **0.7291** | 0.4922 |
| patient-to-patient cos | 0.9999 | **0.7261** | 0.4946 |
| z_biology matrix rank | 16/16 | 16/16 | 16/16 |

The collapse is gone: patient-to-patient cosine holds at 0.726 instead of running to 0.9999, and
cross-modal positives and negatives are now clearly separated (0.9992 vs 0.7291) where before they
were indistinguishable. The loss descends **below chance**.

The fix recovers essentially all of the benefit of deleting decorrelation outright (2.0875 vs 2.0789)
**while keeping the term active**, which is what preserves D1's arm symmetry — `programme_only` and
`programme_free` still differ only in programme supervision.

**G2.6 still does not pass.** The criterion is `biology_contrastive <= 0.10`; we are at 2.0875. A
step-budget/learning-rate sweep (5,000 steps at lr 1e-3; 3,000 at lr 3e-3) is running to establish
whether the remaining gap is optimisation budget or a third defect.

### In plain terms
The model was collapsing because one of its own training rules could be switched off by making every
patient look identical — so it did exactly that. Adding back the rule that forbids patients looking
identical stops the collapse, and the model starts learning to match images to their RNA. It is not
yet learning *well enough* to pass the check, which needs near-perfect recall of 16 patients.

### Meaning for the claim
- The collapse defect is fixed and the fix is verified on real data, not just in a toy.
- The gate threshold was not touched and does not need to be.
- D1 remains blocked pending the sweep; the fix is necessary and demonstrably not yet sufficient.
- Two defects have now been found stacked on this one gate (queue self-cancellation, then collapse).
  A third is possible and the sweep is the cheapest way to find out.

### Files / commits
`v2/losses.py` (hazard warning on `feature_decorrelation`), `v2/training.py` (variance active in both
D1 arms, scale-aware `target_std`, head-scoped application), `v2/tests/test_programme_free.py`
(+2 regression tests), commit `a260dee`.
