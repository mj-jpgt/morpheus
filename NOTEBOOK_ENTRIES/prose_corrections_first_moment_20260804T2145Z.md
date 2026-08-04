## Eight absolute claims withdrawn across five documents: the adjustment discharges the confound from the first moment, a nonlinear probe recovers 3.15×/3.45× chance, and a saturated cell design bounds the residual at 0.6052 → 0.6051

**Logged:** 2026-08-04 21:45 UTC.
**Question:** three documents assert that the confound is *gone* after adjustment. A nonlinear probe
recovers both confounds from the adjusted state. Which sentences are refuted, and what replaces them?
**How obtained:** no new computation. Every number below is quoted from
`NOTEBOOK_ENTRIES/tcga_nonlinear_confound_probe_result_20260804T2100Z.md` and
`NOTEBOOK_ENTRIES/nonlinear_adjustment_channel_result_20260804T2130Z.md`, both read at their committed
state, and each is attributed to the entry and section it came from. Suite re-run on the box at
`~/ws_ind` with threads capped: `test_paper_paths_resolve`, `test_track1_controls`,
`test_calibra_purity`, `test_provenance` — 35 passed.

---

### 1. What was withdrawn, and where

| document | line(s) | withdrawn | replaced by |
|---|---|---|---|
| `v2/research/rebase/nature/PHASE1_RESULT.md` | 41 | "Cancer is gone." | first-moment statement + 3.45× nonlinear reading + the two bounds |
| " | 30 | "survives **full** confound adjustment" | survives the *upper bound* on conditional-mean adjustment (0.6052 → 0.6051) |
| `paper/P1_CALIBRA_DRAFT.md` | §4.2 closing | "no adjusted number in this paper is reading site" | new item **(3)**: first moment, the probe, the two bounds, the anti-prediction finding |
| " | 28 (abstract) | "the adjustment works … zero breaching axes" (unqualified) | "works on the first moment", with the probe and the bound |
| " | 65 (one-para abstract) | "removes what it claims to" | "removes what it claims to **from the first moment**" + the probe + the bound |
| " | 167 (§1.4 contribution 1) | contribution stated without its reach | + "what that certifies is the first moment"; must not be read as "the confound is gone" |
| " | 1428 (§6) | "removes what it claims to" | + "and no further" + the probe + the regenerated null |
| " | 1487 (§6) | "showed the adjustment discharging it" | "discharging it **in the mean**", + what it does not cover, + the anti-prediction finding |
| `paper/P1_FIGURES.md` | 47–48 (F2 claim) | "removes what it claims to remove" | + "from the first moment" + the residual |
| " | 59 (panel d) | bare 0.463 → 0.035 bar | + caption constraint; **new panel (e)** carrying the probe and the bounds |
| `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` | 54, 57 | "**fully discharges** it"; "no adjusted number on this project is reading site" | first-moment + probe + both bounds + **new item 3**, the anti-prediction finding |
| `NOTEBOOK_ENTRIES/t13_adjusted_certificate_and_p6_20260803T0300Z.md` | 1, 27, 34, 89, 91, 106–107 | six absolute sentences | **appended correction block**; history not edited |

The t13 file is an append-only record, so nothing in it was altered: a `## CORRECTION APPENDED` block
at the end tabulates the six withdrawn sentences, states that every *number* in the entry stands
unchanged, and declares itself binding where it contradicts the text above.

### 2. The replacement claim, which is stronger than what it replaces

Written into every corrected document in the same three moves:

1. **The adjustment removes the confound from the first moment.** LDA is a mean-based scorer, so an
   LDA certificate certifies exactly that and no more. The certificate's numbers are untouched:
   0.3633 → 0.0118 at chance 0.0118, a 21–45× drop, zero breaching axes in six states; cancer-type
   0.463 → 0.035 at chance 0.048.
2. **A nonlinear probe still recovers both confounds: 3.15× chance for site, 3.45× for cancer**, netted
   against a null that regenerates the adjustment inside every permutation, *p* at the 1/201 floor,
   three probe families agreeing (k-NN, random forest, RBF-SVM). The mechanism is demonstrated on an
   equal-means synthetic where LDA reads 0.231 and the per-axis maximum 0.244 against chance 0.250
   while k-NN reads 0.554, the forest 0.625 and the SVM 0.658.
3. **The residual is bounded and cannot explain the result.** A **saturated cancer × site cell design**
   — one free parameter per occupied cell — spans every function of the confound labels and therefore
   upper-bounds *any* conditional-mean adjustment, kernel, forest or boosted; it moves the
   `d2_h::wsi_biology` channel by 0.0001, **0.6052 → 0.6051** (retention 0.998). The **confound labels
   alone** reach a channel of 0.1237 (additive) and 0.0903 (saturated), both below the real channel's
   own null median of 0.1483 — **11.2%** and **6.0%** of its excess over that null.

### 3. Recorded for the first time: the certificate's raw row anti-predicts what survives

D1-B `programme_free` (`d1_f`) has the **lowest** raw joint site LDA of the twelve artifacts measured —
0.1071–0.1449, against `d1_p`'s 0.1778–0.3764 and `d2_i_seed43`'s 0.4735 — and the **highest**
post-adjustment nonlinear reading of the twelve, **6.21–7.38× chance for site and 4.98–6.69× for
cancer**. A low raw joint LDA does not predict a low post-adjustment reading; on this cohort it
anti-predicts it. Written into `TRACK1_NEGATIVE_CONTROLS.md` §T1.3 as a new numbered item, into
`P1_CALIBRA_DRAFT.md` §4.2(3) and §6, and into the t13 correction block. Source:
`tcga_nonlinear_confound_probe_result_20260804T2100Z.md` §8 item 2.

### 4. A disagreement between two sources, reported and not resolved

The two entries express the null correction with **different estimators**, and they are not
interchangeable:

| entry | estimator | site | cancer |
|---|---|---:|---:|
| 21:00 §7b | **difference**: observed × chance − regenerated-null median × chance (4.80 − 1.65; 4.67 − 1.22) | **3.15×** | **3.45×** |
| 21:30 §6 | **ratio**: `corrected_multiple` = observed ÷ regenerated-null median | 2.73× (`d2_h` incumbent), range 2.4–3.3× | 3.66× (`d2_h` incumbent), range 3.1–3.7× |

The 3.15/3.45 pair is what the 21:00 entry states it would defend, and it is what the corrected prose
quotes — but every corrected document names the estimator and cites the alternative alongside it, so a
reader cannot mistake one for the other. No document quotes a bare "3.15×" without the qualifier.

A second collision worth naming: **0.231** appears in the 21:00 entry as the equal-means synthetic's
LDA reading and, independently, in the 21:30 entry as the `none` arm's retention of excess. The
corrected prose uses it only in the first sense and always beside its chance rate of 0.250.

### 5. Not edited

* `NOTEBOOK.md:1278` carries this entry's index row with the same absolute verb ("**adjustment fully
  discharges site**"). `NOTEBOOK.md` is out of scope for this pass; flagged here and in the t13
  correction block, not changed.
* `v2/calibra/claim_guards.py` — out of scope, and it contains no string matching any withdrawn claim
  (checked).
* The 21:00 and 21:30 entries themselves, and the other predeclarations and results that already
  *record* the refutation rather than asserting the refuted claim
  (`spatial_claim_replication_result_20260804T1930Z.md`, `p4_certification_end_to_end_20260804T2000Z.md`,
  `paper/P3_P4_PLAN.md:486`) — these are the corrective record and are correct as written.

### Honest constraints

No new measurement was made for this pass; it is a prose reconciliation against two existing entries,
and it inherits every constraint they declare — one cohort, one partition, `wsi_biology` only, the
ceiling bounding functions of the *labels* only, and the corrected null at 200 permutations. The
6.0–11.2% figures are ratios the 21:30 entry states in prose in §5; `labels_only_ceiling` in
`v2/calibra/nonlinear_adjustment.py` returns the components they are formed from but no module
constant holds the ratios themselves, so the corrected documents quote them with their source section
rather than importing them.
