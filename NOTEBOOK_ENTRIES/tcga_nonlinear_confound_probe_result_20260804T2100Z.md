## 2026-08-04 21:00 UTC — TCGA answers the question the spatial run raised: YES. A nonlinear probe recovers both confounds from the adjusted state — tissue source site at 4.5–4.9× chance and cancer type at 3.4–4.9× chance, permutation *p* at the floor — on representations whose certificate reads *at or below* chance

**Logged:** 2026-08-04 21:00 UTC. **Predeclared:**
`NOTEBOOK_ENTRIES/PREDECLARED_tcga_nonlinear_confound_probe_20260804T1840Z.md`, committed `5732c7d`,
**before** any number below existed. **How obtained:** `v2/calibra/nonlinear_confound_probe.py` (new,
14 tests) on the A100 box (150.136.45.194), workspaces `~/ws_probe{,2}` deployed by
`git -c core.autocrlf=false archive HEAD` and verified file-by-file against `git ls-tree -r HEAD`
blob SHA-1 (645 / 646 / 648 files, **0 mismatches**, every time). Thread caps
`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, process parallelism, CPU only (the GPU was on the ALCHEMIST
download; the box carried a co-tenant load of 20–75 on 30 cores throughout, which is why the wall
times below are 2–4× the uncontended cost). Outputs
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/nonlinear_probe/`.

---

### 0. The answer, in one paragraph

**Yes.** After the exact `cancer + pooled TSS` cross-fitted adjustment CALIBRA applies before every
channel number in P1, an out-of-fold k-nearest-neighbour vote still names **which of 85 tissue source
sites** a patient came from **5.3–5.8% of the time against a measured chance rate of 1.18%** (4.5–4.9×),
and still names **which of 21 cancer types** **16.0–23.2% of the time against a measured chance rate of
4.76%** (3.4–4.9×). Every one of those readings has permutation *p* at the **1/201 resolution floor**
against a *global* permutation null. On the same rows, the same folds and the same features, the
certificate's own joint LDA reads **0.0051–0.0098** for site and **0.0382–0.0432** for cancer — **at or
below chance in every case**.

The estimator-family argument therefore transfers from HEST to TCGA. The **magnitude** does not:
HEST's adjusted kNN read 0.729 against a 0.0769 chance (9.5×) — TCGA's reads 4.5–4.9× for site and
3.4–4.9× for cancer. Against the predeclared bands that is **reading C (between 2× and 5×): report the
magnitudes, do not adjudicate** — and every value landed just under the 5× line that would have
triggered reading B. What is *not* band-dependent, and what this entry treats as the finding, is that
the sentences **"the site signal is gone"** and **"Cancer is gone"** are absolute claims and are
refuted at *p* ≤ 0.005 by a classifier that is not a function of class means.

---

### 1. What was run, and what each probe is for

`certify_axes` scores exactly two classifiers — `nearest_class_mean_oof` per axis and
`lda_oof_balanced_accuracy` jointly — and both are functions of the class **means**. The adjustment it
certifies, `cross_fitted_residuals`, is a ridge regression on a one-hot `cancer + pooled TSS` design,
which removes the class **mean vector** by construction. Passing is therefore close to arithmetic. The
three probes added here were chosen so that none of them can be satisfied that way:

* **k-NN vote**, k ∈ {1, 3, 5, 10, 15, 25, 50}, out-of-fold, 5 stratified folds, seed 42. The decision
  rule is the labels of the nearest training rows; no class mean enters at any step. The sweep is not
  a hyperparameter search — **the reading is taken at the maximum over k**, fixed in the
  predeclaration, so that a clean number at one k cannot be selected afterwards and so that a
  representation which only looks clean at large k (where the vote smooths toward a density estimate)
  is caught. Both a plain majority vote and an **inverse-training-frequency weighted** vote are run,
  because the pooled `OTHER` site class is 829 of 2,766 rows (30.0%) and a plain vote can look clean
  purely by collapsing onto it. The plain branch is pinned by a test to reproduce
  `hest_claims.knn_balanced_accuracy_oof` bit-for-bit, so the bulk and spatial readings come from one
  estimator.
* **Random forest**, 300 trees, `max_features="sqrt"`, `class_weight="balanced_subsample"`. Chosen
  because it fails differently from k-NN and shares none of its assumptions: its rule is a set of
  axis-aligned thresholds, it is not a metric method, it is invariant to monotone rescaling of any
  single coordinate, and it can key on a difference in **variance**, in skew, or in an **interaction**
  between two coordinates. All three are invisible to LDA and survive mean-removal untouched.
* **RBF-kernel SVM**, `C=1`, `gamma="scale"`, `class_weight="balanced"`, per-fold standardised. Smooth
  and global where k-NN is local and discrete: the decision function is a weighted sum of Gaussians
  centred on training rows, so it reads the same local geometry without k-NN's failure mode on small
  classes, where a vote tie is resolved by class size.

`v2/tests/test_nonlinear_confound_probe.py` builds classes with **exactly equal means and different
variances** and asserts that on it the certificate's two classifiers sit at chance (LDA 0.231,
per-axis max 0.244, chance 0.250) while the k-NN reads 0.554, the forest 0.625 and the SVM 0.658. That
is the situation the whole question is about, and it is now a test rather than an argument.

**Arms**, all three, every block: `raw`; `adjusted` = `cross_fitted_residuals(state, confound_design(
cancer + pooled TSS))`, seed 42, `n_splits=5`, `alpha=1.0`; and `adjusted_standardised`, the
certificate's own per-axis scaling. A test pins the reconstruction: the joint LDA and per-axis maximum
of `adjusted_standardised` equal `certify_axes(..., residualise=True)` **exactly**, so "the probe reads
the same adjusted state the published numbers come from" is checked, not asserted.

---

### 2. The nulls — and which one applies

Both were computed for every block, by the same imported function, differing only in the strata
argument. `global_permutations` permutes with a single constant stratum and therefore **measures**
chance rather than assuming `1/n_classes`; `within_stratum_permutations` is the certificate's own
within-cancer convention.

**Measured chance matches design chance, so the probe is not capacity-bound.** Global-null medians:
**0.0109–0.0148** for site against a design chance of 0.011765, and **0.0466–0.0492** for cancer
against 0.047619. Global-null p95: 0.0118–0.0226 (site), 0.0519–0.0628 (cancer). Predeclared distrust
item 7 — "if the measured global-null p95 is far above 1/n_classes the probe is too weak here" — is
**discharged**: it is not.

**The within-cancer null is not a chance rate here, and the reason is different from HEST's.**
`nesting_diagnostic` on the test partition, computed before any feature was touched:

| quantity | TCGA test partition | HEST test partition (for contrast) |
|---|---:|---:|
| classes / strata | 85 sites / 21 cancers | 13 slides / 7 oncotree labels |
| strata containing exactly one class | **1 of 21** | 5 of 13 |
| fraction of rows whose label is unpermutable | **0.0101** | 0.385 |
| classes contributing to more than one stratum | **1 of 85** (only pooled `OTHER`) | — |
| `within_stratum_chance` (oracle knows the cancer, then guesses uniformly) | **0.1739** | — |
| design chance `1/n_classes` | 0.011765 | 0.0769 |

So the HEST failure — *a stratum with one class, so the permutation is the identity and the null is
handed the true labels* — is **nearly absent on TCGA** (1.0% of rows, not 38.5%). The live problem is
the other one: **site nests totally inside cancer** (0 of the 84 kept sites contributes patients to two
cancers; only the pooled `OTHER` class spans them), so every within-cancer permuted label is still a
site *of the correct cancer* and the null inherits the entire cancer→site restriction. An oracle that
knew only the cancer would score **0.1739 — 14.8× the design chance**. Quoting that as "chance" would
be wrong by more than an order of magnitude.

**Which applies, as predeclared:** the adjustment under test removes cancer **and** site, so a state
from which both have been regressed out should not be able to name the site at all. **The applicable
bar for the adjusted arms is the global null.** The within-cancer null answers the narrower question
the certificate poses about a *raw* axis ("does it carry site beyond what cancer already explains?"),
and it is reported beside every site row because it decomposes the finding — see §4.

---

### 3. The headline blocks — `wsi_biology`, test partition, n = 2,766

`~/e0_run/d2_v3/d2_v3_s42/artifacts/{d2_h,d2_i}_seed42.npz`, 200 permutations per null, seed 42.
"kNN max" is the maximum over the seven k and the two vote rules, as predeclared. "×" is the multiple
of the design chance rate. *p* is floored at 1/201 = 0.0050.

**Tissue source site — 85 pooled classes, chance 0.011765**

| artifact | arm | joint LDA | per-axis max | **kNN max** | **×** | at | global p95 | global *p* | within-cancer p95 | within-cancer *p* |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| d2_h | raw | 0.3643 | 0.0515 | **0.1511** | **12.84** | prior k50 | 0.0209 | **0.0050** | 0.1100 | 0.0050 |
| d2_h | **adjusted** | **0.0051** | 0.0121 | **0.0554** | **4.71** | prior k3 | 0.0175 | **0.0050** | 0.0490 | 0.0100 |
| d2_h | adjusted, standardised | 0.0052 | 0.0121 | **0.0580** | **4.93** | prior k5 | 0.0188 | **0.0050** | 0.0521 | 0.0199 |
| d2_i | raw | 0.3035 | 0.0432 | **0.1406** | **11.95** | prior k5 | 0.0179 | **0.0050** | 0.0903 | 0.0050 |
| d2_i | **adjusted** | **0.0098** | 0.0121 | **0.0534** | **4.54** | prior k5 | 0.0179 | **0.0050** | 0.0372 | 0.0050 |
| d2_i | adjusted, standardised | 0.0087 | 0.0121 | **0.0576** | **4.90** | prior k3 | 0.0169 | **0.0050** | 0.0386 | 0.0050 |

**Cancer type — 21 test classes, chance 0.047619** (the confounder P1 §4.2 quotes as 0.463 → 0.035)

| artifact | arm | joint LDA | per-axis max | **kNN max** | **×** | at | global p95 | global *p* |
|---|---|---:|---:|---:|---:|---|---:|---:|
| d2_h | raw | 0.6407 | 0.1909 | **0.4289** | **9.01** | prior k25 | 0.0627 | **0.0050** |
| d2_h | **adjusted** | **0.0389** | 0.0464 | **0.2097** | **4.40** | k1 | 0.0541 | **0.0050** |
| d2_h | adjusted, standardised | 0.0399 | 0.0464 | **0.2318** | **4.87** | k1 | 0.0550 | **0.0050** |
| d2_i | raw | 0.5566 | 0.1561 | **0.3927** | **8.25** | prior k50 | 0.0628 | **0.0050** |
| d2_i | **adjusted** | **0.0382** | 0.0442 | **0.1602** | **3.37** | prior k50 | 0.0619 | **0.0050** |
| d2_i | adjusted, standardised | 0.0432 | 0.0442 | **0.1705** | **3.58** | prior k25 | 0.0603 | **0.0050** |

**The k sweep, adjusted arms.** This is the check that decides whether the reading is a k=1 artefact.

| block | k=1 | k=3 | k=5 | k=10 | k=15 | k=25 | k=50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| d2_h cancer, plain | 0.2097 | 0.1884 | 0.1910 | 0.1903 | 0.1906 | 0.1826 | 0.1616 |
| d2_h cancer, prior-corrected | 0.2097 | 0.1961 | 0.1871 | 0.1947 | 0.1820 | 0.1952 | 0.1799 |
| d2_h site, plain | 0.0505 | 0.0450 | 0.0385 | 0.0305 | 0.0266 | 0.0174 | 0.0124 |
| d2_h site, prior-corrected | 0.0505 | 0.0554 | 0.0527 | 0.0332 | 0.0298 | 0.0360 | 0.0347 |

**Cancer is flat across a fifty-fold change in k** (0.16–0.21): it is not a nearest-neighbour
curiosity, it is structure the representation carries at every scale. **Site decays with k under the
plain vote and does not under the prior-corrected one** — which is exactly why predeclared item 6
existed: at k=50 the plain vote reads 0.0124 (1.05× chance, indistinguishable from chance) while the
prior-corrected vote on the same neighbours reads 0.0347 (2.95×, *p* = 0.0050). A run reporting only
the plain vote at k=15 would have concluded "site is gone" from a probe that had collapsed onto the
30% `OTHER` class.

---

### 4. The decomposition the two nulls buy: most of what survives is **cancer**, not site-beyond-cancer

Read the site rows against both nulls at once. Against the **global** null every adjusted value is at
the *p*-floor. Against the **within-cancer** null the adjusted values are only marginally out
(*p* = 0.0050–0.0199, and the observed 0.0534–0.0580 sits just above a within-cancer p95 of
0.0372–0.0521). Meanwhile the **cancer** target reads 0.16–0.23 at *p* = 0.0050.

The honest reading of the pair: **the dominant thing that survives the adjustment to a nonlinear
reader is cancer type, and the surviving site recovery is largely — not entirely — inherited from it.**
Site-beyond-cancer survives too, but weakly. That decomposition is only available because both nulls
were computed; a run quoting the within-cancer null alone would have understated the finding by a
factor of ~4, and a run quoting the global null alone would have attributed to *site* an effect that
is mostly *cancer*.

---

### 5. The other probe families

**Random forest, 300 trees** (observed, `d2_h_seed42` `wsi_biology`, site): **raw 0.0863 (7.3×),
adjusted 0.0222 (1.9×)**. **RBF-kernel SVM** (observed, same block): **raw 0.1957 (16.6×), adjusted
0.0453 (3.85×)**. Measured cost per out-of-fold fit on the contended box: k-NN 0.1 s, SVM 11.5–15.5 s,
forest 126.7 s (raw) / 245.9 s (adjusted).

Two things follow and both were predeclared.

* **Disagreement between families resolves upward** (predeclared distrust item 4). The SVM is the most
  *powerful* reader of the raw state (16.6× chance, above the kNN's 12.8×) and reads 3.85× on the
  adjusted state; the forest is the *weakest* on both (7.3× / 1.9×). Three families spanning metric,
  kernel and axis-aligned decision rules all read the adjusted state above chance, and none of them
  reads it as low as the certificate does.
* **The forest's 1.9× is the one number that could be read as reassuring, and it should not be.** It
  is the family with the lowest power on the *raw* state too (7.3× against 12.8× and 16.6×), so its
  low adjusted reading is a statement about the forest at n = 2,766 with 85 classes, not about the
  representation. Predeclared item 1 — no power, no finding — cuts both ways.

*The permutation nulls for these two families are the expensive part of this run and are reported
separately in §8; the forest's declared 100 permutations were cut to the predeclared floor of 50 for
the reason the predeclaration named in advance (measured wall cost), and the SVM kept its 100.*

---

### 6. A provenance defect found on the way, reported because it was found

**P1 §4.2's joint-LDA numbers do not reproduce from the artifacts everyone believes they came from.**
Three distinct copies of `d2_h_seed42.npz` exist on persistent storage, with three different SHA-256:

| path | mtime | sha256 (first 16) | raw joint LDA | raw per-axis max | adjusted joint LDA | adjusted per-axis max |
|---|---|---|---:|---:|---:|---:|
| `runs/d2_final/artifacts/` | 2026-08-01 20:33 | `4a18b94f1017b85d` | **0.3633** | **0.0532** | **0.0118** | **0.0123** |
| `e0_run/d2_v3/recovered_artifacts/` | 2026-08-03 01:11 | `053490d685bf0dc4` | 0.1782 | 0.0545 | 0.0063 | 0.0120 |
| `e0_run/d2_v3/d2_v3_s42/artifacts/` | 2026-08-03 07:57 | `e81f4496f82c503a` | 0.3785 | 0.0515 | 0.0052 | 0.0121 |
| **published in §4.2** | — | — | **0.3633** | **0.0532** | **0.0118** | **0.0123** |

`confound_certificate.py` has not changed since `942d3c2` (2026-08-02 21:38), i.e. the instrument is
identical; the artifacts are not. **§4.2's row reproduces to four decimal places, on all four
statistics, only from `runs/d2_final/artifacts/`** — while `~/e0_run/d2_v3/*/artifacts/` is the path
the project (and the brief for this run) names. The two differ enough to matter: raw joint LDA 0.3633
against 0.3785, adjusted 0.0118 against 0.0052.

Nothing in §3–§5 depends on which copy is used — the probe was run on both and the conclusion is the
same — but **the paper should name the artifact path and its hash beside the §4.2 table**, and the
`d2_v3` re-export should be recorded as a *different* artifact rather than assumed to be the same one.
This is consistent with the project's own documented finding that training on this stack is not
seed-reproducible; what is new is that a *published* table's provenance path is not the one that
reproduces it.

---

### 7. Every check I said I would run if the result came out favourable — run anyway

The predeclaration listed eight ways a *favourable* result would be untrustworthy. The result is not
favourable, so the symmetric list (§5 of the predeclaration) is the load-bearing one, and both are
discharged here.

1. **Power on the raw arm.** Raw kNN reads 11.95–12.84× chance for site and 8.25–9.01× for cancer, at
   *p* = 0.0050. The probe has power. Reading D is not in play.
2. **Below-chance is a symptom, not a reassurance.** The *certificate's* adjusted numbers are the ones
   sitting below chance (site 0.0051–0.0098 against 0.0118; cancer 0.0382–0.0432 against 0.0476). That
   is the anti-correlation cross-fitted residualisation is known to induce against the variable being
   removed — a statement about the first moment. The kNN, which reads no moment, is *above* its own
   measured null on the same rows.
3. **Max-over-k, fixed in advance.** Reported; and the full sweep is tabled in §3 so the reader can see
   it is not a k=1 artefact for cancer.
4. **Families resolve upward.** Done in §5.
5. **Scale sensitivity.** `adjusted` and `adjusted_standardised` agree in direction and differ by
   ≤ 0.03 absolute (site) and ≤ 0.03 absolute (cancer); the standardised arm is slightly *higher*
   everywhere. Both are reported; neither is quoted alone.
6. **Prior correction.** Reported in both branches. §3 shows the conclusion depends on it at large k
   for site, and does not for cancer.
7. **Capacity.** Discharged in §2 — the measured global-null median equals the design chance.
8. **Patient-key leakage.** `duplicate_patient_ids = 0` in the analysis partition, asserted before any
   probe ran and recorded in every output record. A kNN cannot be naming the site by finding the same
   patient twice.
9. **Cancer standing in for site.** This is the symmetric check, and it *bites*: see §4. The site
   finding is reported as "mostly inherited from surviving cancer" for exactly this reason.
10. **Residual folds versus probe folds.** `cross_fitted_residuals` uses `KFold(seed=42)`; the probe
    uses `_stratified_folds(seed=42)`. These are different partitions, so a probe test row can be a
    residual train row. This is a real non-independence and is recorded rather than waved away. It
    cannot manufacture the effect — the residualisation is a *removal*, and the global permutation null
    is computed on the same already-residualised features, so any structure the residualiser itself
    introduced is inside the null as well as inside the observed.

---

### 8. Breadth, and what did not run

*(filled by the amendment commit)*

---

### 9. Prose that is now wrong, flagged and **not** edited

Per the rules for this run, `NOTEBOOK.md`, the paper drafts and `claim_guards.py` were not touched.
Four things in the current text are unsupportable as written.

1. **`NOTEBOOK_ENTRIES/t13_adjusted_certificate_and_p6_20260803T0300Z.md`, line 27: "This is not
   partial attenuation; the site signal is gone."** Not supported. The supportable sentence is: *the
   site confound is removed from the first moment, and a mean-based certificate therefore certifies; a
   k-NN on the same adjusted rows still recovers it at 4.5–4.9× chance with p at the permutation
   floor.*
2. **Same file, line 34: "Every CALIBRA channel number in this project is measured on adjusted states
   and is therefore **not** reading tissue source site."** Not supported as an absolute. The
   supportable version names the estimator: adjusted states carry no site information *in their class
   means*; they carry site information a nonlinear reader recovers at ~4.7× chance, and cancer
   information it recovers at ~4.4× chance.
3. **`paper/P1_CALIBRA_DRAFT.md` §4.2, closing sentence: "The defect is therefore a property of the raw
   representation, and no adjusted number in this paper is reading site."** Same correction, at the
   prominence of the original claim — this is the paper's most-cited methodological result. The same
   sentence appears at `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md:57`.
4. **`v2/research/rebase/nature/PHASE1_RESULT.md:41`: "cancer-type balanced accuracy from the
   residualised representation drops to 0.035 (chance 0.048) from 0.463 raw. **Cancer is gone.**"** The
   strongest of the four to correct, because the *cancer* number is the larger effect: a k-NN on the
   adjusted representation names the cancer type 16–23% of the time, 3.4–4.9× chance, flat across
   k ∈ [1, 50], at *p* = 0.0050. "Cancer is gone" is false for any reader that is not a class mean.
   §4.2's own summary line (`P1_CALIBRA_DRAFT.md:28`, `:167`, `:1428`) and `paper/P1_FIGURES.md:59`
   carry the same figure and inherit the caveat.

Two additions rather than corrections:

* **The certificate schema should carry the probe family it was issued under.** It already records
  `certified_on = {raw | adjusted}` (P4). "Certified" without "by a mean-based classifier" is the
  ambiguity this whole run exists to remove, and `p4_certify.py`'s `verdict` field inherits it.
* **`confound_certificate.py`'s module docstring argues correctly for the within-cancer null and does
  not say when it stops being a chance rate.** On TCGA it is 0.1739 against a design chance of 0.0118.
  The docstring should carry `nesting_diagnostic`'s two failure modes; the code for that now exists in
  `nonlinear_confound_probe.py` and is tested.

---

### 10. Honest constraints on every number above

* **n = 2,766 patients, 85 site classes** is ~33 patients per class. The probes are not capacity-bound
  (§2) but they are not powerful either; a larger cohort would very likely read higher, not lower.
* **This is one partition of one cohort.** The test partition is the one §4.2 quotes and is the only
  one used for the headline; the 32-class `--partition all` arm was declared and is reported in §8.
* **The k-NN and the SVM are metric probes on a 256-dimensional block.** They read Euclidean structure;
  a confound expressed in a way neither Euclidean geometry nor axis-aligned thresholds can see would
  still be invisible to all three families. "Above chance" is a lower bound on what survives, never an
  upper one.
* **HEST and TCGA are not exchangeable.** HEST spots physically overlap and 13 slides carry the whole
  cohort; TCGA patients do not overlap and 85 sites carry 2,766 of them. That the direction replicates
  is the transferable claim; that the magnitude is 2× smaller is not evidence about either cohort's
  quality, and the spatial 1c caveat stays attached to the spatial numbers.

---

### 11. Suite status

At the commit this entry lands on, on the deployed workspace: **`morpheus/v2/tests` is 421 passed, 27
errors** in `~/venv` (the 27 being `test_p2_figures.py`, which needs matplotlib, absent from that venv
by policy), and **28 passed** for `test_p2_figures.py` in `~/axis_venv` (matplotlib 3.10.9). The 14 new
tests in `v2/tests/test_nonlinear_confound_probe.py` are inside the 421; the pre-existing baseline was
408 + the other agents' additions.

### 12. Files, commits, outputs

* `v2/calibra/nonlinear_confound_probe.py`, `v2/tests/test_nonlinear_confound_probe.py` — commits
  `5cfccf4`, `a170d57`, `1992f1a` on `research/rebase-vision`.
* Predeclaration `5732c7d`.
* Results JSON: `p1_evidence/nonlinear_probe/{anchor_knn,anchor_svm,d2final_wsi,breadth_knn}.json`
  under `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/`.
* Logs: `p1_evidence/nonlinear_probe/logs/*.log`.
