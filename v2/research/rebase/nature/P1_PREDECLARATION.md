# P1 evidence — PREDECLARATION (written before any Track 2 / dilution run)

Written 2026-08-03 UTC on the Lambda box, workspace `~/ws_p1`, before a single sweep or
dilution cell was executed in this session. Every prior run output was lost when
`/home/ubuntu` was wiped on instance stop; the code survived in the repo snapshot, the
numbers did not. So all numbers reported in TRACK2_INDUCED_CORRELATION.md and
DILUTION_LOWER_BOUND.md are produced after this file, not before it.

## Provenance caveat that is NOT covered by this predeclaration

`v2/calibra/induced_correlation_sweep.py` already contained, on disk, four candidate
forms (P0, P1, P2, P3). Its own docstring states that **P3 was written after seeing P1
and P2 fail on the rank ladder at n=2,530 / seed 42**. P3 is therefore *post hoc on the
rank axis* and is reported as such. It is out of sample on every other n and every other
seed, and that is the only sense in which it is tested here. P0/P1/P2 are genuinely
a-priori. Nothing below rehabilitates P3.

## B — induced correlation

Quantity: `|r_induced|` = median over 40 draws of `corr(X_res u, A_res v)` at spike
level 0, where `A` is constructed exactly orthogonal to `X u` before residualisation.
Anchor: 0.067-0.140 for the 99-column cancer+TSS design at n=2,530.

### Predictions, in the order they will be graded

* **P0 (the plan text)** — `|r_induced| ~ k/n`; log-log exponents `(b_k, b_n) = (+1, -1)`.
  Grows with design rank, vanishes with n.
* **P1** — `0.6745 * kappa / sqrt(k_eff)`, `kappa = R_s R_a / sqrt((1-R_s^2)(1-R_a^2))`.
* **P2** — `0.6745 * R_s R_a / sqrt(k_eff)`.
* **P3 (post hoc on rank; see caveat)** — `0.6745 * R_x R_y / sqrt(k_eff_shared)`.
* **P4 (mine, new, a-priori)** — **the effect is a bias, not a sampling fluctuation.**
  `|r_induced|` converges to a non-zero constant as n grows. Concretely: at the anchor
  design, `|r_induced|(n=6,427) >= 0.6 * |r_induced|(n=2,530)`.
  **Falsifier: a ratio below 0.6, with the ratio tracking sqrt(2530/6427)=0.63 or lower
  across all real designs, would mean the effect is a finite-sample artefact and P0 is
  closer to right than the mechanism I am asserting.**
* **P5 (mine, a-priori, the decisive falsifier)** — the effect is *structural*: it needs a
  design that predicts BOTH modalities. A matched-rank design that predicts neither must
  induce nothing. Concretely: for `gaussian_k99` and `permuted_cancer_tss_pool10` at
  n=2,530, `|r_induced| <= 0.025`, i.e. at most ~2x the pure sampling scale
  `0.6745/sqrt(n) = 0.013`, and at least 3x smaller than the real design at the same rank
  and n. **Falsifier: a matched-rank structureless design reproducing 0.067-0.140 would
  mean the phenomenon is degrees-of-freedom bookkeeping (already fully covered by the
  classical N-R effective-sample-size result) and our magnitude claim would carry no
  content beyond it.**
* **P6 (estimator, a-priori)** — varying `n_splits` over {2, 5, 10, 20} and Ridge `alpha`
  over {0.01, 1.0, 100.0} moves the anchor `|r_induced|` by less than 25% relative.
  **Falsifier: >2x movement means we are measuring our residualiser, not residualisation.**

### What would sink Track 2 entirely
P5 failing. If a structureless matched-rank design gives the same magnitude, the only
content left is the textbook N-R result and Track 2 becomes a replication note.

## C — Winkler (settled BEFORE B was run; recorded here for order-of-work)
Verdict: **NO magnitude.** Winkler et al. 2020 NeuroImage 220:117065 reports error rates
(pcer/FWER, %), power (%) and p-value calibration, never a correlation-unit magnitude for
residualisation-induced dependence between independent blocks. Detail and quotes in
TRACK2_INDUCED_CORRELATION.md.

## D — dilution lower bound

Arm: `foreign_tumour` — same-cancer, different-patient tumour patches from the same
H-Optimus store, drawn donor-slide-first, nested across levels
d in {0.0, 0.10, 0.20, 0.40, 0.60}. Representation: `raw_hoptimus_meanstd`
(concat(mean, std), no fitted parameters). Readout: CALIBRA on the same 99-column
cancer+TSS design, same partition, same seed.

### Predictions
* **D1** — `heldout_top_cca` (the channel) declines **monotonically** in d.
  Falsifier: any non-monotone step larger than its own bootstrap CI.
* **D2** — the decline is approximately proportional to the surviving own-tumour
  fraction: `channel(d) / channel(0) ~= (1-d)` to within +/-0.15 absolute over
  d <= 0.40. Falsifier: a ratio outside [(1-d)-0.15, (1-d)+0.15].
* **D3** — the unpaired `detection_floor` is non-decreasing in d.
* **D4** — at d = 0.60 the channel retains **more than 0** but **less than 55%** of
  its d = 0 value; i.e. contamination is costly but not annihilating.

### The "lower bound" claim, and why it is an assumption, not a result
The claim is: foreign same-cancer tumour patches are the *most benign* possible
non-informative contaminant (same stain, same preparation, same cancer, same store), so
the measured cost is a floor on the cost of any other contaminant. **This is an
assumption I cannot test without the normal-tissue arm, and there is a mechanism arguing
the other way**: normal tissue is off-manifold in a way that is similar across ALL
patients, so it adds a near-constant offset to a mean-pooled bag, which damages
between-patient variation less than adding a *patient-specific random* tumour shift does.
If that mechanism dominates, `foreign_tumour` is an upper bound, not a lower one.
**Declared in advance: the write-up will state the direction as assumed-and-unverified,
and will not assert "lower bound" without this caveat attached.**

## Rules binding this session
* Every headline number carries a permutation null at >= 1,000 permutations and a CI.
  `perm p = 1/(n+1)` is reported as a resolution floor, never as "p < 0.05".
* Must-FAIL controls go through `GateLedger.add`; baseline comparisons through
  `GateLedger.observe`.
* Any must-FAIL control that passes is reported as an instrument defect.
