## 2026-08-03 19:30 UTC — G2.6 passes but `programme_free` is collapsing in real training: effective rank 1.76 against `programme_only`'s 7.38

**Logged:** 2026-08-03 19:30 UTC. **How obtained:** `~/ws_d1/geom_probe.py` on the A100
(`150.136.45.194`), loading the live D1 checkpoints at epoch 21–23 and measuring `z_biology` on
**282 held-out test patients**. Training metrics from `~/e0_run/d1_v1/*/train_metrics.jsonl`.
**The run has not been stopped** — see "what was not done" below.

### Technical

Measured on held-out test patients, mid-training:

| arm | epoch | centred eff-rank | hard rank @1e-3 | RNA–RNA cos | feat-std |
|---|---|---|---|---|---|
| `programme_only` seed 42 | 21 | **7.38** | 104 | 0.317 | 0.0044 |
| `programme_only` seed 43 | 21 | **7.35** | 115 | 0.317 | 0.0048 |
| `programme_free` seed 42 | 23 | **1.76** | **9** | **0.977** | 0.0137 |

The two `programme_only` seeds agree to two decimal places, so the estimate is stable rather than
noisy. `programme_free`'s biology representation is **effectively 1.76-dimensional across 282
patients**, with a hard rank of 9, and its RNA-view biology states sit at mutual cosine 0.977 — very
nearly one vector.

**Read the raw cosine carefully.** Uncentred WSI–WSI cosine is ~0.99 for *both* arms. That is *not*
by itself pathological and must not be quoted as one: `z_biology` is normalised but never centred, so
a large batch-common direction is expected (it is ~81% of the squared norm at initialisation — the
finding that drove today's G2.6 work). The meaningful quantity is the geometry *after* centring,
which is what a linear channel sees and what CALIBRA residualises to. There the arms separate
decisively: 7.38 versus 1.76.

**The signature is rank-1 collapse, not shrinkage.** `programme_free` has the *higher* per-feature
std (0.0137 vs 0.0044) and the *lower* rank. That is exactly the `z_i = m + a_i·u` family — plenty of
variance along a single direction — which is the failure mode established earlier today, and the
precise reason a per-dimension variance floor cannot prevent it.

**Onset coincides with the end of warmup.** From the training trajectory, `programme_free`'s
per-feature std rises to 0.041 by epoch 6 and then falls monotonically to 0.015–0.018, while its
contrastive term drifts *upward* from 6.15 toward its chance value of ln(4310) = 8.369. Warmup ends
at epoch 4, which is when `decorrelation` (0.04) and `variance` (0.01) switch on. `programme_only`
over the same span moves the other way — per-feature std 0.033 → 0.049, close to the isotropic 0.0625.

**Why this is most likely a defect and not the answer to D1's question.** Today's controlled
measurements on the same architecture and the same contrastive objective:

* with `decorrelation` excluded, the memorisation check reaches contrastive 0.012–0.057, retrieval
  16/16, patient cosine 0.0597 and effective rank 5.81;
* with `decorrelation` at 0.04 it sits at 2.61 against a chance of 2.7726, and at *every* weight from
  0.001 to 4.0 it collapses the representation within 25 steps while switching itself off
  (20.74 → 0.00);
* `feature_decorrelation` has total collapse as its global minimum, which is already pinned by
  `test_decorrelation_is_minimised_by_the_collapse_it_claims_to_prevent`.

Additionally, the RNA biology view — the most collapsed at 0.977 — carries **no** anti-collapse term
at all: `decorrelation` is applied only to `out_wsi["z_biology"]`, and the variance floor only to the
fused `output["z_biology"]`. The one view with no regulariser is the one that collapsed hardest.

**Consequence for D1 as configured.** If this holds to epoch 40, D1 compares a 7.4-dimensional
representation against a 1.8-dimensional one, and `programme_only` will win the molecular channel for
a reason that has nothing to do with programme supervision. The pre-registered prediction is
`programme_free >= programme_only`, and the manifest's instruction for `programme_only` winning is to
escalate rather than reframe — but that instruction anticipated a scientific outcome, not a collapsed
arm. **This would be a third outcome the gate did not enumerate, and quoting it as evidence about
supervision would be wrong.**

**This also bears on D2, and that is not my call to make.** D2's arms are both `programme_only`-family
and ran with the same `--decorrelation-weight 0.04`. `programme_only` is not collapsed here — 7.38 is
low but far from degenerate — so this is not a claim that D2 is invalid. It is a flag that D2's
measured channel sits on a representation whose centred effective rank is ~7 out of 282, which is
worth knowing when interpreting a channel difference.

### What was NOT done, deliberately

**The run was not stopped and the objective was not changed mid-flight.** Three reasons:

1. I could be wrong about the endpoint. Contrastive training can recover, and stopping at epoch 22 on
   a projection would destroy a possibly-valid run to save three GPU-hours on an otherwise idle box.
2. The completed run *is* the evidence. Audit check A5 measures effective rank per arm per seed
   precisely so this is quantified at epoch 40 rather than argued from mid-training.
3. Changing `decorrelation` changes the experiment's design and both arms' regularisation. That is a
   scientific decision, and the honest move is to escalate with evidence rather than to quietly
   reconfigure an experiment while it runs.

If the collapse holds at epoch 40, the recommendation is: re-run D1 with `--decorrelation-weight 0`
in **both** arms (symmetric, so the contrast still isolates programme supervision), and treat this
run as the control demonstrating why. That would require changing
`test_both_d1_arms_pair_decorrelation_with_a_variance_floor`, which asserts `decorrelation > 0` in
both arms — a test encoding a theory that today's measurements have falsified.

### In plain terms

The health check now passes, and the training that follows it is still going wrong — for the same
underlying reason, in a place the health check does not look. The check deliberately switches off two
penalty terms because they interfere with measuring whether the model can memorise; real training
leaves them on, and one of them is the thing that crushes the representation.

Concretely: the arm trained without curated supervision has ended up describing every patient with
essentially one number, where the supervised arm uses about seven. If that holds, comparing them
tells you which arm collapsed, not which kind of supervision works better.

I have not stopped the run or changed the recipe. Both are decisions with scientific consequences,
the evidence will be stronger at epoch 40 than at epoch 22, and the machine is otherwise idle.

### Meaning for the claim

D1 is at risk of producing an uninterpretable answer for a reason that is now well characterised and
fixable. Nothing should be quoted from it until the epoch-40 effective ranks are read. The
biology-head-trained-without-programme-supervision arm that `PHASE1B` says exists nowhere on disk may
still not exist in usable form after this run.

### Files / commits

- `~/ws_d1/geom_probe.py`, `~/e0_run/d1_v1/*/train_metrics.jsonl`, `~/e0_run/d1_v1/*/last.pt`
- `v2/losses.py` `feature_decorrelation` docstring; `v2/training.py` decorrelation/variance placement
- `v2/tests/test_programme_free.py::test_decorrelation_is_minimised_by_the_collapse_it_claims_to_prevent`
- Prior entries: `g26_rank_collapse_diagnosis`, `g26_term_isolation`, `g26_passes`
