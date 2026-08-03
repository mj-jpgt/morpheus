## 2026-08-03 01:52 UTC — T1.3 must-FAIL control 1 DID NOT FAIL. Every D2 state predicts tissue source site far above its own permutation null; no axis on either artifact can currently be certified

**Logged:** 2026-08-03 01:52 UTC. **How obtained:** `python -m morpheus.v2.calibra.confound_certificate --artifacts d2_h_seed42.npz d2_i_seed42.npz --partition test --n-permutations 1000 --n-boot 200 --n-boot-axes 8 --n-jobs 6`, Lambda box `~/ws_p1`, 2,766 held-out patients, 85 pooled TSS classes (`min_site_count=10`), within-cancer label permutation null, seed 42.

### Technical

This is the control that says: **no axis we would certify may be able to predict which hospital the
slide came from.** It is condition 4 of P4's five-point certification rule, verbatim. The pass
criterion was fixed in advance: per-axis out-of-fold balanced accuracy for pooled TSS must not exceed
the 95th percentile of a ≥1,000-draw within-cancer label-permutation null, and the axis bootstrap CI
must include the chance rate.

**Result: it did not fail. Every state, on both artifacts, is NOT CERTIFIED and fails the joint test.**
Chance rate = 1/85 = 0.0118. Permutation resolution = 1/1001 = 0.000999.

| artifact | state | per-axis max bal. acc. | per-axis median | median null p95 | axes breaching (of 256) | **joint LDA bal. acc.** | joint null p95 | joint perm p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| d2_h | wsi_biology | 0.0532 | 0.0334 | 0.0392 | **17 (6.6%)** | **0.3633** | 0.1539 | 0.000999 |
| d2_h | full_biology | 0.0548 | 0.0322 | 0.0359 | **60 (23.4%)** | 0.2630 | 0.1742 | 0.000999 |
| d2_h | rna_biology | 0.0506 | 0.0320 | 0.0360 | **58 (22.7%)** | 0.2563 | 0.1744 | 0.000999 |
| d2_i | wsi_biology | 0.0511 | 0.0298 | 0.0345 | **43 (16.8%)** | 0.2348 | 0.1224 | 0.000999 |
| d2_i | full_biology | 0.0551 | 0.0334 | 0.0367 | **61 (23.8%)** | 0.2689 | 0.1856 | 0.000999 |
| d2_i | rna_biology | 0.0495 | 0.0326 | 0.0362 | **48 (18.8%)** | 0.2744 | 0.1861 | 0.000999 |

Every joint permutation p is at the 1/1001 resolution floor: **not one of 1,000 within-cancer label
permutations reached the observed joint accuracy for any state.**

**The shape of the leak matters as much as its existence.** Individual axes are weak site predictors —
the best single axis anywhere reaches 0.055, only 4.6× chance, and the median axis (0.030–0.033) sits
*below* its own permutation null p95. But the **joint** linear discriminant over all 256 axes reaches
0.235–0.363, i.e. **20–31× chance** and roughly 2× the joint null p95. Site is not concentrated in a
few nameable axes; it is smeared across the representation in a combination that no per-axis screen
would catch. That is precisely the failure mode the joint test was added for, and it is the reason a
per-axis-only certificate would have issued a clean bill of health here.

**What this closes.** `v2/calibra/e0_basis_transfer.py:923` currently records G3.5 as
`unavailable_no_site_labels`. It is no longer unavailable and it is no longer clean: TSS labels are
derivable from the patient barcode and, once derived, the gate does not pass.

**What is still pending, and it is the decisive follow-up.** The adjusted arm — the identical
certificate run on `cross_fitted_residuals(state, cancer + pooled TSS)`, i.e. after the exact
adjustment CALIBRA applies before it measures any channel — is still running
(`p1_evidence/track1/certificate_adjusted/`). The two arms answer different questions and the
difference between them decides how much of this is a P1 problem versus a P3/P4 problem:

* raw arm (this entry) → *may this axis be exposed to a user?* Answer: **no**.
* adjusted arm (pending) → *does our own adjustment actually discharge the site confound, or only
  appear to?* If the adjusted arm certifies, every CALIBRA channel number in the project remains
  interpretable and the defect is confined to the raw representation. If the adjusted arm also fails,
  the adjustment is not removing what it claims to remove and **every adjusted number on the project
  is reading site**.

The certificate's own self-test covers this direction: `test_certificate_residualisation_removes_a_pure_site_code`
asserts the adjusted accuracy drops by >0.2 on a planted pure site code, so the machinery is known to
be able to show a discharge when one occurs.

### In plain terms

We built the check that is supposed to prove our axes are not secretly reading "which hospital sent
this slide". It came back the wrong way. No single axis is much of a hospital detector — the best one
is right about 5.5% of the time where guessing would be 1.2%. But if you let a simple classifier use
all 256 axes together, it identifies the source hospital about a third of the time out of 85
possibilities. Out of a thousand random relabellings, none came close. So the hospital signal is
genuinely there; it is just spread thin enough that looking at axes one at a time would have missed
it entirely.

This is a defect and it is being recorded as one, not explained away. It does not by itself invalidate
the channel measurements, because those are made *after* the site adjustment — but whether the
adjustment works is exactly the run that has not finished yet.

### Meaning for the claim

* **P4: no axis of any D2 state may be exposed.** Condition 4 of the five-point certification rule
  fails for all six state/artifact combinations. Any `legible_axis` claim built on these artifacts is
  currently inadmissible.
* **P3:** every per-axis legibility claim on D2 inherits the same block.
* **P2:** the "effective rank went up" comparison can now be met with "your representation is partly a
  site code" and, on the raw states, that objection is correct. P2 must quote the *adjusted* arm.
* **P1 (this paper):** this is the battery working. A must-FAIL control that fails to fail is the most
  valuable output a negative-control battery can produce, and it is worth more to P1 than a clean pass
  would have been — provided it is reported as a defect rather than tuned away. It also demonstrates
  the structural point the battery exists to make: the joint test catches what per-axis screening
  cannot, so a certificate without it is decoration.
* **A design lesson for the certificate spec itself:** the T1.3 criterion as written in the plan is
  per-axis only. Had it been implemented literally, `wsi_biology` on d2_h would have shown 17/256 axes
  breaching and a per-axis maximum of 0.053 — arguably dismissible — while the real leak is 31× chance
  in the joint direction. The joint test must be a required field of the certificate schema, not an
  optional extra.

### Files / commits

`v2/calibra/confound_certificate.py` (inherited, audited, unchanged).
Results `p1_evidence/track1/certificate_raw/{confound_certificate.json,task_rows.csv}`;
adjusted arm pending at `p1_evidence/track1/certificate_adjusted/` — both under
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/`.
