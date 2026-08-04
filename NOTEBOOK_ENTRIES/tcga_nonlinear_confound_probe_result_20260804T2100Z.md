## 2026-08-04 21:00 UTC — TCGA answers the question the spatial run raised: YES. On the exact artifact P1 §4.2 quotes, a k-NN recovers tissue source site at **4.80× chance** and cancer type at **4.67× chance**, and an RBF-SVM recovers cancer at **6.10× chance**, all at the permutation floor — from states whose certificate reads **0.0118** and **0.0458**, i.e. at or below chance

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

### 0. The answer

**Yes.** On `d2_h_seed42` `wsi_biology`, test partition, n = 2,766 — the exact block P1 §4.2 quotes, on
the artifact that reproduces §4.2's row to four decimal places — the **adjusted** state whose joint LDA
reads **0.0118** (chance 0.0118) is read by an out-of-fold k-NN at **0.0565 balanced accuracy for
tissue source site, 4.80× chance**, permutation *p* = 0.0050 at the 1/201 floor against a **measured**
global null whose p95 is 0.0156. The same state's cancer-type joint LDA reads **0.0458** (chance
0.0476, i.e. below chance) while a k-NN reads **0.2224, 4.67× chance**, *p* = 0.0050, and an RBF-kernel
SVM reads **0.2904, 6.10× chance**, *p* = 0.0099 at its own 1/101 floor.

Both confounders survive, in both artifacts, on both artifact generations, and in **all three** probe
families. On that same §4.2 state a **random forest** — axis-aligned thresholds, not a metric method —
reads site at **0.0335, 2.85× chance** and cancer at **0.2696, 5.66× chance**, both at 2.7–5.0× their
own global null p95 and at the 1/51 floor. The certificate cannot see any of it because it scores class
means and the adjustment removes class means.

**Against the predeclared bands** (§3 of the predeclaration, stated in terms of the k-NN): every
adjusted k-NN reading is between 2× and 5× chance, which is **reading C — report the magnitudes, do not
adjudicate**. Two things must be said beside that, and both were fixed in advance:

* The predeclaration's distrust item 4 says **disagreement between families resolves upward**. On the
  cancer target **two of the three families are over the 5× line** — the SVM at **6.10×** and the
  forest at **5.66×** — and the k-NN, the family the bands were written in terms of, is at 4.67×. The
  letter of the band is C; the substance is on the reading-B side for every family except the one the
  bands name.
* Reading C is about **how much** survives. It is not about **whether** anything does. The sentences
  *"the site signal is gone"* and *"Cancer is gone"* are absolute, and they are refuted at
  *p* ≤ 0.005 by three classifiers that are not functions of class means. That correction does not
  depend on which band the magnitude falls in.

**It survives the strictest null anyone has proposed.** A sharper null published mid-run by another
agent (`regenerated_adjustment_null`, `efee0f8`) regenerates the adjustment inside every permutation,
so each draw carries an adjustment artefact of its own size. It shows the artefact is real — it lifts
the reference from 1.0× to **1.65×** chance for site and **1.22×** for cancer — and that the finding
outlives it: observed values sit at 2.1× and 3.2× that null's p95, *p* at the floor, leaving a netted
**3.15× chance for site and 3.45× for cancer**. §7b, which also withdraws a defence I made and got
wrong.

The estimator-family argument therefore transfers from HEST to TCGA; the **magnitude** does not. HEST's
adjusted k-NN read 0.729 against a 0.0769 chance (9.5×, 73% of spots). TCGA's reads 4.3–4.8× for site
(5% of patients) and 3.6–4.9× for cancer (17–22% of patients).

---

### 1. What was run, and why each probe

`certify_axes` scores exactly two classifiers — `nearest_class_mean_oof` per axis and
`lda_oof_balanced_accuracy` jointly — and both are functions of the class **means**. The adjustment it
certifies, `cross_fitted_residuals`, is a ridge regression on a one-hot `cancer + pooled TSS` design,
which removes the class **mean vector** by construction. Passing is close to arithmetic. The three
probes were chosen so none of them can be satisfied that way:

* **k-NN vote**, k ∈ {1, 3, 5, 10, 15, 25, 50}, out-of-fold, 5 stratified folds, seed 42. The rule is
  the labels of the nearest training rows; no class mean enters at any step. The sweep is not a
  hyperparameter search — **the reading is the maximum over k**, fixed in the predeclaration, so a
  clean number at one k cannot be chosen afterwards and a representation that only looks clean at large
  k (where the vote smooths toward a density estimate) is caught. Both a plain majority vote and an
  **inverse-training-frequency weighted** vote are run, because the pooled `OTHER` site class is 829 of
  2,766 rows (30.0%) and a plain vote can look clean purely by collapsing onto it. The plain branch is
  pinned by a test to reproduce `hest_claims.knn_balanced_accuracy_oof` bit-for-bit, so the bulk and
  spatial readings come from one estimator.
* **Random forest**, 300 trees, `max_features="sqrt"`, `class_weight="balanced_subsample"`. Chosen
  because it fails differently from k-NN and shares none of its assumptions: its rule is a set of
  axis-aligned thresholds, it is not a metric method, it is invariant to monotone rescaling of any
  single coordinate, and it can key on a difference in **variance**, in skew, or in an **interaction**
  between two coordinates. All three are invisible to LDA and survive mean-removal untouched.
* **RBF-kernel SVM**, `C=1`, `gamma="scale"`, `class_weight="balanced"`, per-fold standardised. Smooth
  and global where k-NN is local and discrete: the decision function is a weighted sum of Gaussians
  centred on training rows, so it reads the same local geometry without k-NN's failure mode on small
  classes, where a vote tie is resolved by class size. Because it standardises inside each fold it is
  scale-invariant by construction, which is why its `adjusted` and `adjusted_standardised` readings are
  identical below — a useful consistency check rather than a duplicate.

`v2/tests/test_nonlinear_confound_probe.py` builds classes with **exactly equal means and different
variances** and asserts that on it the certificate's two classifiers sit at chance (LDA 0.231,
per-axis max 0.244, chance 0.250) while the k-NN reads 0.554, the forest 0.625 and the SVM 0.658. That
is the situation the whole question is about, and it is a test rather than an argument.

**Arms**, all three, every block: `raw`; `adjusted` = `cross_fitted_residuals(state, confound_design(
cancer + pooled TSS))`, seed 42, `n_splits=5`, `alpha=1.0`; and `adjusted_standardised`, the
certificate's own per-axis scaling. A test pins the reconstruction: the joint LDA and per-axis maximum
of `adjusted_standardised` equal `certify_axes(..., residualise=True)` **exactly**. On the real data
that identity is visible in the tables below — `adjusted_standardised` reproduces §4.2's published
0.0118 / 0.0123 (d2_h) and 0.0052 / 0.0104 (d2_i) exactly, so the probe demonstrably reads the state
the published numbers come from.

---

### 2. The nulls — and which one applies

Both were computed for every block by the same imported function, differing only in the strata
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
| `within_stratum_chance` (oracle knows the cancer, then guesses uniformly among its sites) | **0.1739** | — |
| design chance `1/n_classes` | 0.011765 | 0.0769 |

The HEST failure — *a stratum with one class, so the permutation is the identity and the null is handed
the true labels* — is **nearly absent on TCGA** (1.0% of rows, not 38.5%). The live problem is the
other one: **site nests totally inside cancer** (0 of the 84 kept sites contributes patients to two
cancers; only the pooled `OTHER` class spans them), so every within-cancer permuted label is still a
site *of the correct cancer* and the null inherits the whole cancer→site restriction. An oracle knowing
only the cancer would score **0.1739 — 14.8× the design chance**. Quoting that as "chance" would be
wrong by more than an order of magnitude.

**Which applies, as predeclared:** the adjustment under test removes cancer **and** site, so a state
from which both have been regressed out should not be able to name the site at all. **The applicable
bar for the adjusted arms is the global null.** The within-cancer null answers the narrower question
the certificate poses about a *raw* axis, and it is reported beside every site row because it
decomposes the finding — §4.

---

### 3. The headline: `runs/d2_final/artifacts/`, the artifacts that reproduce P1 §4.2 exactly

`wsi_biology`, partition `test`, n = 2,766, 200 permutations per null, seed 42. "kNN max" is the
maximum over the seven k and the two vote rules, as predeclared. "×" is the multiple of the design
chance rate. *p* is floored at 1/201 = 0.0050.

**Tissue source site — 85 pooled classes, chance 0.011765**

| artifact | arm | joint LDA | per-axis max | **kNN max** | **×** | at | global p95 | global *p* | within-cancer p95 | within-cancer *p* |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| d2_h | raw | 0.3494 | **0.0532** | **0.1560** | **13.26** | prior k15 | 0.0217 | **0.0050** | 0.1091 | 0.0050 |
| d2_h | adjusted | 0.0147 | 0.0123 | **0.0511** | **4.34** | k3 | 0.0154 | **0.0050** | 0.0399 | 0.0050 |
| d2_h | **adjusted, standardised** | **0.0118** | **0.0123** | **0.0565** | **4.80** | k3 | 0.0156 | **0.0050** | 0.0428 | 0.0050 |
| d2_i | raw | 0.2237 | **0.0511** | **0.1261** | **10.72** | prior k25 | 0.0223 | **0.0050** | 0.0950 | 0.0050 |
| d2_i | adjusted | 0.0061 | 0.0104 | **0.0507** | **4.31** | k1 | 0.0155 | **0.0050** | 0.0335 | 0.0050 |
| d2_i | **adjusted, standardised** | **0.0052** | **0.0104** | **0.0548** | **4.65** | prior k3 | 0.0164 | **0.0050** | 0.0419 | 0.0050 |

The bolded `adjusted, standardised` rows are the **published §4.2 numbers, reproduced**: joint LDA
0.0118 / per-axis max 0.0123 for d2_h and 0.0052 / 0.0104 for d2_i, against §4.2's table. (The raw rows
show per-axis max 0.0532 and 0.0511, also §4.2's published values; the raw *joint* LDA differs from
§4.2's 0.3633 / 0.2348 only because this arm is unstandardised — `certify_axes` run directly on the
same file returns 0.3633, see §6.)

**Cancer type — 21 test classes, chance 0.047619** (the confounder §4.2 quotes as 0.463 → 0.035)

| artifact | arm | joint LDA | per-axis max | **kNN max** | **×** | at | global p95 | global *p* |
|---|---|---:|---:|---:|---:|---|---:|---:|
| d2_h | raw | 0.6238 | 0.1956 | **0.4466** | **9.38** | prior k25 | 0.0620 | **0.0050** |
| d2_h | adjusted | 0.0442 | 0.0480 | **0.2126** | **4.46** | prior k5 | 0.0582 | **0.0050** |
| d2_h | adjusted, standardised | 0.0458 | 0.0480 | **0.2224** | **4.67** | prior k5 | 0.0595 | **0.0050** |
| d2_i | raw | 0.4584 | 0.1599 | **0.3649** | **7.66** | prior k50 | 0.0624 | **0.0050** |
| d2_i | adjusted | 0.0357 | 0.0466 | **0.1738** | **3.65** | prior k15 | 0.0609 | **0.0050** |
| d2_i | adjusted, standardised | 0.0391 | 0.0466 | **0.1852** | **3.89** | prior k10 | 0.0611 | **0.0050** |

### 3b. Replication on the second artifact generation, `~/e0_run/d2_v3/d2_v3_s42/artifacts/`

Same states, same protocol, different export of the same seed (see §6). Site, adjusted: kNN max
**0.0554 (4.71×)** d2_h and **0.0534 (4.54×)** d2_i; standardised, **0.0580 (4.93×)** and
**0.0576 (4.90×)**. Cancer, adjusted: **0.2097 (4.40×)** and **0.1602 (3.37×)**; standardised
**0.2318 (4.87×)** and **0.1705 (3.58×)**. Raw: site **0.1511 (12.84×)** / **0.1406 (11.95×)**, cancer
**0.4289 (9.01×)** / **0.3927 (8.25×)**. Every global *p* = 0.0050. The conclusion does not depend on
which export is used.

### 3c. The k sweep — the check that decides whether this is a k = 1 artefact

`d2_v3_s42` `d2_h`, adjusted arms, balanced accuracy:

| block | k=1 | k=3 | k=5 | k=10 | k=15 | k=25 | k=50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| cancer, plain vote | 0.2097 | 0.1884 | 0.1910 | 0.1903 | 0.1906 | 0.1826 | 0.1616 |
| cancer, prior-corrected | 0.2097 | 0.1961 | 0.1871 | 0.1947 | 0.1820 | 0.1952 | 0.1799 |
| site, plain vote | 0.0505 | 0.0450 | 0.0385 | 0.0305 | 0.0266 | 0.0174 | 0.0124 |
| site, prior-corrected | 0.0505 | 0.0554 | 0.0527 | 0.0332 | 0.0298 | 0.0360 | 0.0347 |

**Cancer is flat across a fifty-fold change in k** (0.16–0.21): not a nearest-neighbour curiosity,
structure the representation carries at every scale. **Site decays with k under the plain vote and does
not under the prior-corrected one** — which is exactly why predeclared item 6 existed: at k = 50 the
plain vote reads 0.0124 (1.05× chance, indistinguishable from chance) while the prior-corrected vote on
the same neighbours reads 0.0347 (2.95×, *p* = 0.0050). A run reporting only the plain vote at k = 15
would have concluded "site is gone" from a probe that had collapsed onto the 30% `OTHER` class.

---

### 4. What the two nulls buy: most of what survives is **cancer**, not site-beyond-cancer

Read the site rows against both nulls at once. Against the **global** null every adjusted value is at
the *p*-floor and 3–4× the null p95. Against the **within-cancer** null the adjusted values are also at
the floor but much closer to it (observed 0.0507–0.0565 against a within-cancer p95 of 0.0335–0.0428).
Meanwhile the **cancer** target reads 0.17–0.22 at *p* = 0.0050 with a k-NN and 0.29 with an SVM.

The honest reading of the pair: **the dominant thing that survives the adjustment to a nonlinear reader
is cancer type; the surviving site recovery is real but is substantially inherited from it.** That
decomposition is only available because both nulls were computed. A run quoting the within-cancer null
alone would have understated the site effect ~3×; a run quoting the global null alone would have
attributed to *site* an effect that is largely *cancer*.

---

### 5. The other two families

**RBF-kernel SVM**, 100 permutations, global and within-cancer nulls, `d2_v3_s42` `wsi_biology`,
*p* floored at 1/101 = 0.0099:

| artifact | target | arm | SVM balanced accuracy | × chance | global p95 | global *p* |
|---|---|---|---:|---:|---:|---:|
| d2_h | site | raw | 0.1957 | **16.63** | 0.0212 | **0.0099** |
| d2_h | site | adjusted (= standardised) | 0.0453 | **3.85** | 0.0188 | **0.0099** |
| d2_i | site | raw | 0.1807 | **15.36** | 0.0189 | **0.0099** |
| d2_i | site | adjusted (= standardised) | 0.0427 | **3.63** | 0.0184 | **0.0099** |
| d2_h | cancer | raw | 0.5076 | **10.66** | 0.0626 | **0.0099** |
| d2_h | cancer | **adjusted (= standardised)** | **0.2904** | **6.10** | 0.0610 | **0.0099** |

**The SVM is the most powerful reader of both the raw and the adjusted state for cancer**, and its
adjusted cancer reading of 6.10× chance is above the predeclared reading-B line. It is also
scale-invariant by construction — identical on `adjusted` and `adjusted_standardised` — so this is not
a standardisation artefact.

**Random forest**, 300 trees, 50 permutations (the predeclared cost floor), global null, on the exact
§4.2 block `runs/d2_final/artifacts/d2_h_seed42.npz` `wsi_biology`, *p* floored at 1/51 = 0.0196:

| target | arm | joint LDA | forest balanced accuracy | × chance | global null p95 | global *p* |
|---|---|---:|---:|---:|---:|---:|
| site | raw | 0.3494 | 0.0755 | **6.42** | 0.0133 | **0.0196** |
| site | **adjusted, standardised** (§4.2's published 0.0118 state) | **0.0118** | **0.0335** | **2.85** | 0.0125 | **0.0196** |
| cancer | raw | 0.6238 | 0.4537 | **9.53** | 0.0519 | **0.0196** |
| cancer | **adjusted, standardised** | **0.0458** | **0.2696** | **5.66** | 0.0536 | **0.0196** |

**The forest is the third family and it agrees — and on cancer it is the second family over the 5×
line.** On the state the certificate scores at 0.0458 against a chance rate of 0.0476, a classifier
built from axis-aligned thresholds — not a metric method, cannot be satisfied by mean-removal — names
the cancer type **27.0% of the time, 5.66× chance, 5.0× its own permutation null p95**. On the site
target it reads 2.85× at 2.7× its null p95. Both at the *p*-floor.

The forest is the *weakest* family on site (2.85× against the k-NN's 4.80× and the SVM's 3.85×) and it
is also the weakest on the site raw arm (6.42× against 13.26× and 16.63×), so that lower site number is
about forest power at 85 classes with ~33 patients each, not about the representation — predeclared
item 1 cutting in the other direction. At 21 classes it has the samples to work with and reads high.

Measured cost per out-of-fold fit on the contended box: k-NN 0.1 s, SVM 11.5–15.5 s, forest 126.7 s
(raw) / 245.9 s (adjusted). That cost is why the forest's declared 100 permutations were cut to the
predeclared floor of **50**, for the reason the predeclaration named in advance, and why the forest
null was run on a reduced block set. No permutation count was ever increased after seeing a result.

---

### 6. A provenance defect found on the way, reported because it was found

**P1 §4.2's numbers come from an artifact path the project does not name.** Three distinct copies of
`d2_h_seed42.npz` exist on persistent storage, with three different SHA-256, and `certify_axes` — whose
file has not changed since `942d3c2` (2026-08-02 21:38), so the instrument is identical — returns three
different answers:

| path | mtime | sha256 (first 16) | raw joint LDA | raw per-axis max | adjusted joint LDA | adjusted per-axis max |
|---|---|---|---:|---:|---:|---:|
| `runs/d2_final/artifacts/` | 2026-08-01 20:33 | `4a18b94f1017b85d` | **0.3633** | **0.0532** | **0.0118** | **0.0123** |
| `e0_run/d2_v3/recovered_artifacts/` | 2026-08-03 01:11 | `053490d685bf0dc4` | 0.1782 | 0.0545 | 0.0063 | 0.0120 |
| `e0_run/d2_v3/d2_v3_s42/artifacts/` | 2026-08-03 07:57 | `e81f4496f82c503a` | 0.3785 | 0.0515 | 0.0052 | 0.0121 |
| **published in §4.2** | — | — | **0.3633** | **0.0532** | **0.0118** | **0.0123** |

**§4.2's row reproduces to four decimal places, on all four statistics, only from
`runs/d2_final/artifacts/`** — while `~/e0_run/d2_v3/*/artifacts/` is the path the project (and the
brief for this run) names. The two differ enough to matter: raw joint LDA 0.3633 against 0.3785,
adjusted 0.0118 against 0.0052.

Nothing in §3–§5 depends on which copy is used — both were run and the conclusion is identical — but
**the paper should name the artifact path and its hash beside the §4.2 table**, and the `d2_v3`
re-export should be recorded as a *different* artifact rather than assumed to be the same one. This is
consistent with the project's own documented finding that training on this stack is not
seed-reproducible; what is new is that a *published* table's provenance path is not the one that
reproduces it.

---

### 7. Every predeclared check, run and reported

The predeclaration listed eight ways a *favourable* result would be untrustworthy and three ways an
*unfavourable* one would be. The result is unfavourable, so the second list is load-bearing; both are
discharged.

1. **Power on the raw arm.** Raw k-NN reads 10.72–13.26× chance for site and 7.66–9.38× for cancer, at
   *p* = 0.0050; the raw SVM reads 15.36–16.63× for site. The probe has power. Reading D is not in play.
2. **Below-chance is a symptom, not a reassurance.** The *certificate's* adjusted numbers are the ones
   at or below chance (site 0.0052–0.0147 against 0.0118; cancer 0.0357–0.0458 against 0.0476). That is
   the anti-correlation cross-fitted residualisation is known to induce against the variable being
   removed — a statement about the first moment. The probes, which read no moment, sit above their own
   measured nulls on the same rows and folds.
3. **Max-over-k, fixed in advance.** Reported, with the whole sweep in §3c so the reader can see it is
   not a k = 1 artefact for cancer.
4. **Families resolve upward.** The SVM's 6.10× for adjusted cancer is the highest reading and is
   reported as the finding, not averaged away.
5. **Scale sensitivity.** `adjusted` and `adjusted_standardised` agree in direction everywhere and
   differ by ≤ 0.006 absolute (site) and ≤ 0.012 absolute (cancer); the standardised arm is slightly
   *higher* every time. The SVM is identical across the two by construction. Both are reported.
6. **Prior correction.** Reported in both branches. §3c shows the site conclusion depends on it at
   large k, and the cancer conclusion does not.
7. **Capacity.** Discharged in §2 — the measured global-null median equals the design chance.
8. **Patient-key leakage.** `duplicate_patient_ids = 0` in the analysis partition, asserted before any
   probe ran and recorded in every output record. A k-NN cannot be naming the site by finding the same
   patient twice.
9. **Cancer standing in for site** (the symmetric check, and it bites): see §4. The site finding is
   reported as substantially inherited from surviving cancer for exactly this reason.
10. **Residual folds versus probe folds.** `cross_fitted_residuals` uses `KFold(seed=42)`; the probe
    uses `_stratified_folds(seed=42)`. Different partitions, so a probe test row can be a residual train
    row. **My first defence of this was wrong and §7b replaces it.**

---

### 7b. The objection to item 10, and the null that settles it

While this run was in flight another agent (`efee0f8`) identified a real asymmetry in the null used
above, and it lands directly on item 10:

> `calibration.permutation_null` — the *channel's* null — re-residualises `y` on every permutation, so
> a correlation the shared residualisation *induces* is regenerated inside the null.
> `probe_state` permutes the **labels of an already-adjusted matrix**. A structure the adjustment tied
> to the *true* rows is therefore broken in the null and intact in the observed, and is scored as
> surviving confound rather than as chance.

That mechanism is concrete: cross-fitted residualisation displaces each (class × residual-fold) group
by that fold's estimation error, an offset shared by every patient in the group. My §7 item 10 claimed
"the residualisation is a removal, so it cannot manufacture the effect". **That claim was too strong
and I am withdrawing it rather than defending it.**

`nonlinear_adjustment.regenerated_adjustment_null` is the null that closes the gap: it permutes the
**rows of the block before the adjustment**, so every draw carries an adjustment artefact of its own
size while labels, design, class structure and probe folds stay fixed, and it takes the max over the
same k grid and both vote rules in the null as in the observed. Run on the exact §4.2 block
(`runs/d2_final/artifacts/d2_h_seed42.npz`, `wsi_biology`, test, `adjusted_standardised` — the state
whose certificate reads 0.0118), 200 permutations, seed 42, *p* floored at 1/201:

| target | chance | observed | × chance | **regenerated-null median** | × chance | regenerated-null p95 | *p* | **excess over the null median, in × chance** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| tissue source site (85) | 0.011765 | 0.0565 | 4.80 | **0.0194** | 1.65 | 0.0266 | **0.0050** | **3.15** |
| cancer type (21) | 0.047619 | 0.2224 | 4.67 | **0.0581** | 1.22 | 0.0699 | **0.0050** | **3.45** |

**The objection is real and the finding survives it.** The adjustment artefact is measurable — it lifts
the reference from 1.0× chance (the label-permutation null, median 0.0113 site / 0.0470 cancer) to
**1.65× and 1.22×** — so **part of the raw multiple-of-chance figure is an artefact of the adjustment,
and any future run must quote the regenerated null, not the label-permutation null.** But the observed
values sit at **2.1× and 3.2× the regenerated null's p95**, with *p* at the resolution floor in both
targets. The honest headline figures, net of the adjustment artefact, are **3.15× chance for site and
3.45× chance for cancer**, and they are the numbers I would defend.

This is also the reason the §0 verdict does not move: the predeclared bands were stated on
multiple-of-chance, where the netted figures are still comfortably inside band C, and the *existence*
claim — that "the site signal is gone" and "Cancer is gone" are false for a non-mean reader — is
established against the strictest null available, not the most convenient one.

---

### 8. Breadth: ten more artifacts, including all six D1-B

Predeclared as "no thresholds attached, reported for breadth": the remaining `d2_v3` artifacts (seeds
43, 44) and the six D1-B artifacts `~/e0_run/d1_v2/artifacts/d1_{p,f}_seed4{2,3,4}.npz`, state
`wsi_biology`, partition `test`, arms `raw` and `adjusted`, **k ∈ {5, 15} only** (plain and
prior-corrected), 200 permutations, both nulls. Because the sweep is truncated these are **lower
bounds** on max-over-k: at the anchor the maximum often sat at k = 1 or k = 3.

**Adjusted arms, k-NN max over the truncated grid, as a multiple of chance.** Every global *p* below
is 0.0050, the floor.

| artifact | site × chance | site within-cancer *p* | cancer × chance |
|---|---:|---:|---:|
| d2_h_seed43 | 4.42 | 0.0746 | 4.56 |
| d2_h_seed44 | 4.32 | 0.0348 | 4.13 |
| d2_i_seed43 | 3.23 | 0.1493 | 3.27 |
| d2_i_seed44 | **5.20** | 0.0050 | 4.28 |
| d1_p_seed42 | 2.69 | 0.5970 | 3.47 |
| d1_p_seed43 | 3.52 | 0.0597 | 3.55 |
| d1_p_seed44 | 4.23 | 0.0846 | 4.11 |
| **d1_f_seed42** | **6.21** | 0.0050 | **5.82** |
| **d1_f_seed43** | **6.67** | 0.0100 | 4.98 |
| **d1_f_seed44** | **7.38** | 0.0050 | **6.69** |

Raw arms for the same blocks read 6.33–13.18× (site) and 6.72–9.71× (cancer), so the probe has power
everywhere.

Three things in this table matter.

1. **Four of these twenty adjusted readings exceed 5× chance**, up to **7.38×** — on a truncated k grid
   and therefore understated. These blocks carry no predeclared threshold (the bands were scoped to the
   anchor), so they do not change the band verdict; they do show the anchor is not the worst case.
2. **D1-B `programme_free` (`d1_f`) is the worst of the twelve artifacts on both confounders**, and it
   is the artifact whose *raw* certificate looks best: its raw joint LDA for site is 0.1071–0.1449
   against `d1_p`'s 0.1778–0.3764 and `d2_i_seed43`'s 0.4735. **A low raw joint LDA does not predict a
   low post-adjustment nonlinear reading — on this cohort it anti-predicts it.** That is a direct
   argument against using the certificate's joint row as a proxy for "how confounded is this
   representation".
3. **The within-cancer *p* column shows the §4 decomposition holding across artifacts**: on five of the
   ten blocks the adjusted site reading is *not* significant against the within-cancer null
   (*p* = 0.06–0.60) while every one is at the floor against the global null. Where site survives
   beyond cancer it does so weakly; what survives strongly is cancer.

**The declared 32-class arm.** `--partition all`, cancer target, 32 classes, chance 0.03125,
`runs/d2_final/artifacts/`, 100 permutations, global null (*p* floor 0.0099), k ∈ {1, 5, 15, 50}:

| artifact | arm | joint LDA | kNN max | × chance | global p95 | global *p* |
|---|---|---:|---:|---:|---:|---:|
| d2_h | raw | 0.6129 | 0.4145 | 13.26 | 0.0404 | 0.0099 |
| d2_h | adjusted | 0.0412 | 0.1865 | **5.97** | 0.0381 | 0.0099 |
| d2_h | adjusted, standardised | 0.0396 | 0.1997 | **6.39** | 0.0371 | 0.0099 |
| d2_i | raw | 0.4426 | 0.3556 | 11.38 | 0.0405 | 0.0099 |
| d2_i | adjusted | 0.0278 | 0.1543 | **4.94** | 0.0392 | 0.0099 |
| d2_i | adjusted, standardised | 0.0281 | 0.1592 | **5.09** | 0.0370 | 0.0099 |

**Caveat stated with the number:** `--partition all` includes the 3,118 training patients, so this arm
is in-distribution and is not a held-out measurement. It is reported because the 32-class figure was
declared in advance and because the project's framing quotes 32 classes; it is **not** used for any
verdict. The held-out 21-class arm in §3 is.

### 8b. What did not run, and why

* **The random forest's permutation null was cut**, exactly as the predeclaration allowed and for the
  reason it named in advance: 126.7 s (raw) / 245.9 s (adjusted) per out-of-fold fit at n = 2,766 with
  85 classes, against 0.1 s for the k-NN and 11.5–15.5 s for the SVM. Its declared 100 permutations
  became the predeclared floor of **50**, and its block set was reduced from the twelve-block anchor
  grid to `d2_h_seed42` `wsi_biology` on `runs/d2_final/`, `raw` and `adjusted_standardised`. **The
  site rows completed and are in §5.** The two **cancer**-target forest blocks were still running when
  this entry was written and are not reported; the forest's cancer reading on TCGA is therefore **not
  measured here**. The corresponding k-NN and SVM cancer rows are, and they are the ones the verdict
  rests on.
* **The other two states** (`full_biology`, `rna_biology`) were run only through `certify_axes`'
  published numbers, not through the probe: the probe grid was spent on `wsi_biology` (the image-only
  channel, which is what §4.2 headlines) across twelve artifacts rather than on three states of two.
  Nothing about the estimator argument is state-specific, but the measurement for those two states on
  TCGA is **not made here** and should not be assumed from `wsi_biology`.
* **`min_site_count`** was left at the project default of 10 throughout; no sensitivity sweep on the
  pooling threshold was run.

---

### 9. Prose that is now wrong, flagged and **not** edited

Per the rules for this run, `NOTEBOOK.md`, the paper drafts and `claim_guards.py` were not touched.
Four things in the current text are unsupportable as written.

1. **`NOTEBOOK_ENTRIES/t13_adjusted_certificate_and_p6_20260803T0300Z.md`, line 27: "This is not
   partial attenuation; the site signal is gone."** Not supported. The supportable sentence is: *the
   site confound is removed from the first moment, and a mean-based certificate therefore certifies; a
   k-NN on the same adjusted rows still recovers it at 4.3–4.9× chance with p at the permutation
   floor.*
2. **Same file, line 34: "Every CALIBRA channel number in this project is measured on adjusted states
   and is therefore **not** reading tissue source site."** Not supported as an absolute. The
   supportable version names the estimator: adjusted states carry no site information *in their class
   means*; they carry site information a nonlinear reader recovers at ~4.8× chance and cancer
   information it recovers at 4.7× (k-NN) to 6.1× (SVM).
3. **`paper/P1_CALIBRA_DRAFT.md` §4.2, closing sentence: "The defect is therefore a property of the raw
   representation, and no adjusted number in this paper is reading site."** Same correction, at the
   prominence of the original claim — this is the paper's most-cited methodological result. The same
   sentence appears at `v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md:57`.
4. **`v2/research/rebase/nature/PHASE1_RESULT.md:41`: "cancer-type balanced accuracy from the
   residualised representation drops to 0.035 (chance 0.048) from 0.463 raw. **Cancer is gone.**"** The
   strongest of the four to correct, because the *cancer* number is the larger effect: on the adjusted
   representation a k-NN names the cancer type 17–22% of the time (3.7–4.7× chance, flat across
   k ∈ [1, 50]) and an RBF-SVM 29% of the time (6.10× chance), both at their permutation floor.
   "Cancer is gone" is false for any reader that is not a class mean. §4.2's own summary line
   (`P1_CALIBRA_DRAFT.md:28`, `:167`, `:1428`) and `paper/P1_FIGURES.md:59` carry the same figure and
   inherit the caveat.

Two additions rather than corrections:

* **The certificate schema should carry the probe family the certificate was issued under.** It already
  records `certified_on = {raw | adjusted}` (P4). "Certified" without "by a mean-based classifier" is
  the ambiguity this run exists to remove, and `p4_certify.py`'s `verdict` field inherits it.
* **`confound_certificate.py`'s module docstring argues correctly for the within-cancer null and does
  not say when it stops being a chance rate.** On TCGA it is 0.1739 against a design chance of 0.0118.
  The docstring should carry `nesting_diagnostic`'s two failure modes; the code for that now exists in
  `nonlinear_confound_probe.py` and is tested.

---

### 10. Honest constraints on every number above

* **n = 2,766 patients over 85 site classes** is ~33 patients per class. The probes are not
  capacity-bound (§2) but they are not powerful either; a larger cohort would likely read higher, not
  lower.
* **One partition of one cohort.** The test partition is the one §4.2 quotes and is the only one used
  for the headline.
* **The k-NN and the SVM are metric probes on a 256-dimensional block.** They read Euclidean structure;
  a confound expressed in a way neither Euclidean geometry nor axis-aligned thresholds can see would
  still be invisible to all three families. "Above chance" is a lower bound on what survives, never an
  upper one.
* **HEST and TCGA are not exchangeable.** HEST spots physically overlap and 13 slides carry the whole
  cohort; TCGA patients do not overlap and 85 sites carry 2,766 of them. That the direction replicates
  is the transferable claim; that the magnitude is ~2× smaller is not evidence about either cohort's
  quality, and the spatial 1c caveat stays attached to the spatial numbers.

---

### 11. Suite status

On the deployed workspace `~/ws_probe3` at HEAD, `morpheus/v2/tests` runs **467 passed, 1 failed, 27
errors** in `~/venv`; the 27 errors are `test_p2_figures.py` needing matplotlib, absent from that venv
by policy, and that file is **28 passed** in `~/axis_venv` (matplotlib 3.10.9). My 14 tests in
`v2/tests/test_nonlinear_confound_probe.py` pass (verified in isolation as well as in the suite).

**The one failure is not mine and is flagged rather than touched:**
`test_inductive_adjustment.py::test_one_row_at_a_time_equals_the_whole_block`, added at `e071d6c`
today by another agent. It asserts `np.allclose(a, b, atol=0, rtol=0)` — bit-exact equality — on two
arrays that print identically to eight decimal places, i.e. it is a last-bit floating-point difference
between the batched and the one-row-at-a-time path. It was failing before my last commit and is
untouched by anything in this run. Owner's call whether the right fix is a tolerance or a genuine
associativity bug in the operator; I have not edited the file.

The tree has grown under other agents throughout this session (645 → 669 tracked files), which is why
the passing count differs from the 408 baseline quoted in earlier entries.

### 12. Files, commits, outputs

* `v2/calibra/nonlinear_confound_probe.py`, `v2/tests/test_nonlinear_confound_probe.py` — commits
  `5cfccf4`, `a170d57`, `1992f1a` on `research/rebase-vision`.
* Predeclaration `5732c7d`.
* Results JSON under
  `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/nonlinear_probe/`:
  `d2final_wsi.json` (headline), `anchor_knn.json`, `anchor_svm.json`, `allpart_cancer.json`,
  `breadth_knn.json`, `d2final_forest.json`, `regenerated_null.json` (§7b). Logs in `logs/`.
* §7b uses `v2/calibra/nonlinear_adjustment.regenerated_adjustment_null`, another agent's instrument
  (`efee0f8`), imported rather than reimplemented; the driver is
  `p1_evidence/nonlinear_probe/logs/regen.log`.
