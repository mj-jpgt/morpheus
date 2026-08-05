## 2026-08-05 01:10 UTC — The inductive channel result was **not** a property of one split, but the *precision* of the two published figures was: across **12** discovery/exposure partitions `d2_h::wsi_biology` retention runs **0.9872–1.0524** (median 1.0102) and `d2_i` runs **0.9221–1.0794** (median 1.0165) — so `d2_i` breaks the predeclared ≤ 0.10 spread bar at **0.157**, the published 0.9966 / 0.9710 are the **second-lowest of twelve** on both artifacts, and the global-null caveat of **0.757** is a lone outlier whose next-worst partition reads 0.9226

**Logged:** 2026-08-05 01:10 UTC. **Predeclared in**
`NOTEBOOK_ENTRIES/PREDECLARED_inductive_channel_split_stability_20260805T0015Z.md`, committed
(`3155f5e`) before the `--split-seed` wiring was written and before any number below existed.

**How obtained.** Workspace `~/ws_p1sp/morpheus` on the A100 (`150.136.45.194`), built from
`git -c core.autocrlf=false archive HEAD` at commit `2b0ad53`. **341/341** tracked
`v2|tests|src|configs × .py|.json|.yaml` files verified against an md5 manifest generated from that
commit: **340 identical, 1 reported differing and traced** — `v2/research/rebase/p2/figures/data/
ws_amp/out/P2_CENTRING_VERDICT.json` is stored in git *with CRLF*, so the workspace's bytes match the
blob exactly (`dc34a840…`) and the "mismatch" is my manifest's LF-normalisation, the same artefact
the 23:45 entry recorded for its `.pdf`/`.png` files. **0 code files differ.** CPU only,
`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1`, `--n-jobs 5–6` on a 30-core box carrying a co-tenant load
of 19–22 throughout; the GPU sat at 99% for other agents and was never touched. Python 3.10.12,
numpy 2.2.6, scikit-learn 1.7.2, scipy 1.15.3, pandas 2.3.3.

Driver `v2/research/rebase/nature/p1_evidence/inductive_channel.py`, the same one, **plus one
default-preserving flag**. Operator `v2/calibra/inductive_adjustment.ConfoundAdjustmentOperator`,
split `p4_certify.exposure_split`, state `p4_certify.prepare_state`, statistics
`nonlinear_adjustment.{channel_under_adjustment, labels_only_ceiling, retention_of_excess,
adjuster_agreement, cross_fitted_r2, _load_block}` — **all imported unchanged. No second operator,
no second split function, no statistic defined inline.** `claim_guards.py`, `claim_evidence.json`,
other agents' `PREDECLARED_*` files and `paper/P1_CALIBRA_DRAFT.md` are **not edited**.

**26 runs**: 12 partitions × 2 artifacts at the published discovery fraction 0.5 (the seed-42 one
being the reproduction gate, and the only run carrying the n = 2,766 arm), plus a deliberate coverage
gradient at f = 0.3 / 0.7 on `d2_h`. Artifacts `runs/d2_final/artifacts/d2_{h,i}_seed42.npz`, state `wsi_biology`, `test`
partition, targets `data/frozen_rna_targets.npz`, the same **90** non-control targets, 16 components,
2,000-permutation within-cancer pairing null, `min_site_count` 10.

---

## Bad news first

**1. `d2_i`'s retention is partition-dependent by more than the predeclared bar allows. N2 fires.**
§4 of the predeclaration voided the point-quotation of retention if `max − min` over the partitions
exceeded 0.10. On `d2_i` it is **0.1573** (0.9221 at split seed 31337, 1.0794 at 555) — and it was
already 0.1573 on the **predeclared eight**, before any extension. `d2_h` clears the bar at
**0.0652** over twelve (0.0311 over the predeclared eight). **So "retention 0.9710" may not be
written as a number. The supportable form is a range with the number of partitions stated.**

**2. The two published figures are the second-lowest of twelve, on both artifacts.** `d2_h`'s 0.9966
is rank 2/12 (only split seed 9's 0.9872 is lower) and `d2_i`'s 0.9710 is rank 2/12 (only 31337's
0.9221 is lower). The direction is the *opposite* of the failure mode this run was hunting — the
seed-42 partition was an unlucky draw, not a lucky one, so N3 does not fire and the published claim
is if anything understated. **It is still a defect: a figure quoted as a point estimate turns out to
be near the bottom of its own distribution, and nobody knew the distribution existed.** In 9 of 12
`d2_h` partitions and 8 of 12 `d2_i` partitions the inductive arm's excess **exceeds** its matched
in-sample control's, i.e. retention > 1.

**3. The `0.757` global-null caveat now in P1 §4.2's provenance footnote is a single-partition
artefact, and a dramatic one.** Global-pairing-null retention on `d2_h` across twelve partitions:
**0.7572, 0.9226, 0.9316, 0.9355, 0.9363, 0.9617, 0.9669, 0.9858, 0.9866, 0.9882, 0.9953, 1.0129** —
median **0.9643**. The seed-42 value is not merely the minimum, it sits **0.165 below the next-worst
partition**, while the remaining eleven span 0.090. On `d2_i` the same figure runs 0.8386–1.0322,
median 0.9115. The paper currently offers 0.757 as what "a reader who prefers the global convention
should read". On eleven of twelve partitions that reader would read 0.92–1.01.

**4. The labels-only ceiling's absolute share is even less stable than the 23:45 entry showed, and
the instability is *worse in the transductive arm*.** Share of the channel's excess, `d2_h`, twelve
partitions, identical encoding and identical n throughout: additive design **transductive
+0.0223 → +0.3777, a 16.9× spread**, against **inductive +0.1792 → +0.4148, 2.3×**. Saturated cell
design: transductive 0.0243 → 0.1920 (7.9×), inductive 0.0213 → 0.4346. `d2_i` reproduces both
patterns (15.0× / 2.3×). **The published 6.0% / 11.2% are not merely encoding-dependent, as the
23:45 entry established — at fixed encoding and fixed n they move by an order of magnitude across
partitions of the same cohort.**

**5. None of the above changes the answer to the question this run was commissioned to ask.** Every
one of 24 partition-arms clears the permutation floor, no partition falls below 0.92, the mechanism
sweep says coverage costs ~0.02 of retention per 14 points of coverage, and the labels-only ceiling
on the encoding the adjustment spans is **insignificant in 12 of 12 partitions on both artifacts**.
**The channel survives inductive adjustment on every partition tested.**

---

## 1. The reproduction gate, exact

Nothing below is comparable to the 23:45 entry without this, and it was run first.

| quantity | 23:45 entry | this run | |
|---|---:|---:|:---:|
| `d2_h` `transductive_full` S1 / null / excess / S2 / rank | 0.6052 / 0.1483 / 0.4569 / 0.5841 / 22.50 | **identical to every digit** | ✔ |
| `d2_h` `none_exposure` S1 / null / excess | 0.8234 / 0.7262 / 0.0972 | **identical** | ✔ |
| `d2_h` `transductive_exposure` S1 / null / excess | 0.6173 / 0.2080 / 0.4092 | **identical** | ✔ |
| `d2_h` `inductive_exposure` S1 / null / excess | 0.6145 / 0.2067 / 0.4078 | **identical** | ✔ |
| `d2_h` retention | 0.9966 | **0.9966** | ✔ |
| `d2_i` full / exposure / inductive S1, retention | 0.4703 / 0.5020 / 0.4898, 0.9710 | **identical**, **0.9710** | ✔ |
| `d2_h` ceilings, tr/ind × additive/saturated/frozen | 35.8 / 27.1 / 19.2 / 18.9 / 6.4 / −19.3 % | **35.84 / 27.09 / 19.20 / 18.90 / 6.45 / −19.31 %** | ✔ |
| operator digest, exposure rows, adjusted blocks vs `prepare_state` | `2060a635…d9c2ce`, bit-identical | **`2060a635…d9c2ce`, bit-identical ×3** | ✔ |
| discovery-fold site coverage | 606 / 1,382 = 43.8% | **606 / 1,382 = 0.43849** | ✔ |

**R1 discharged.** The n = 2,766 arm was run once, in this gate; every other run carries
`--skip-full` because that arm does not depend on the partition.

## 2. The partitions are twelve different partitions — §5.1 of the predeclaration

`split_overlap_probe.py` (imports `exposure_split` and `_load_block`, defines nothing): over the
predeclared eight, pairwise **Jaccard of the exposure patient sets 0.3168 – 0.3450, median 0.3327** —
the ~1/3 two independent halves of a cohort must have — **8 distinct partitions**, every fold
1,384 / 1,382 with **21 cancers on both sides** in every split, **0 patients in both folds** asserted
in the driver on every run. The predeclaration voided the run above 0.9. Nothing near it.

Each run also records `exposure_patient_digest`; the twelve differ, and only the seed-42 run
reproduces P4's operator digest. On every other partition the `prepare_state` bit-identity
comparison is **recorded as skipped with its reason**, not reported as a False.

## 3. The spread — the primary result, reported before the individual numbers

`wsi_biology`, `test`, n = 1,382 exposure rows in every row of both tables, retention =
`nonlinear_adjustment.retention_of_excess(inductive_arm, matched_transductive_control)` on the
project's within-cancer pairing null.

| | `d2_h` predeclared 8 | `d2_h` all 12 | `d2_i` predeclared 8 | `d2_i` all 12 |
|---|---:|---:|---:|---:|
| median retention | 1.0109 | **1.0102** | 1.0065 | **1.0165** |
| min | 0.9966 | **0.9872** | 0.9221 | **0.9221** |
| max | 1.0276 | **1.0524** | 1.0794 | **1.0794** |
| **max − min** | 0.0311 | **0.0652** | **0.1573** | **0.1573** |
| s.d. | 0.0107 | 0.0168 | 0.0482 | 0.0423 |
| retention > 1 | 6/8 | 9/12 | 5/8 | 8/12 |
| seed 42's rank (1 = lowest) | 1/8 | **2/12** | 1/8 | **2/12** |
| every arm at the permutation floor | 8/8 | **12/12** | 8/8 | **12/12** |

**Why twelve and not eight, stated plainly.** The predeclaration committed to eight partitions
({42, 7, 11, 23, 101, 555, 2718, 31337}). After seeing `d2_i` breach the ≤ 0.10 bar on those eight I
added four more ({4242, 9, 77, 1000}) **to characterise the range, not to change a verdict** — and an
extension can only ever widen a min–max range, never narrow it, so it cannot flatter the result. It
did widen `d2_h` (0.0311 → 0.0652, seed 9 and seed 77 becoming the new extremes) and did **not**
widen `d2_i` (both extremes were already in the predeclared eight). Both the predeclared-8 and the
all-12 figures are in the table and the verdict is graded on the eight.

### `d2_h::wsi_biology`, all twelve partitions

| split seed | **retention** | retention (global null) | site coverage | frequent sites | operator design | ind. S1 | ind. null | ind. excess | tr. S1 | tr. null | tr. excess | ind. S2 | tr. S2 | p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 9 | **0.9872** | 0.9363 | 0.4457 | 33 | 57 | 0.6022 | 0.2091 | 0.3931 | 0.6073 | 0.2092 | 0.3982 | 0.5694 | 0.5771 | 0.0005 |
| **42** (published) | **0.9966** | **0.7572** | 0.4385 | 31 | 55 | 0.6145 | 0.2067 | 0.4078 | 0.6173 | 0.2080 | 0.4092 | 0.5701 | 0.5805 | 0.0005 |
| 31337 | 0.9967 | 0.9617 | 0.4870 | 39 | 63 | 0.5695 | 0.2055 | 0.3640 | 0.5752 | 0.2100 | 0.3652 | 0.5133 | 0.5409 | 0.0005 |
| 11 | 1.0024 | 0.9226 | 0.4747 | 35 | 59 | 0.6102 | 0.2101 | 0.4001 | 0.6112 | 0.2120 | 0.3992 | 0.5557 | 0.5406 | 0.0005 |
| 1000 | 1.0090 | 0.9953 | 0.4609 | 36 | 60 | 0.6223 | 0.2022 | 0.4201 | 0.6230 | 0.2067 | 0.4163 | 0.4400 | 0.5328 | 0.0005 |
| 4242 | 1.0098 | 0.9858 | 0.4711 | 38 | 62 | 0.6012 | 0.2073 | 0.3939 | 0.6001 | 0.2101 | 0.3900 | 0.5618 | 0.5569 | 0.0005 |
| 7 | 1.0106 | 0.9866 | 0.4356 | 32 | 56 | 0.6296 | 0.2039 | 0.4257 | 0.6294 | 0.2082 | 0.4212 | 0.6011 | 0.6161 | 0.0005 |
| 23 | 1.0112 | 0.9355 | 0.4595 | 33 | 57 | 0.5783 | 0.2049 | 0.3734 | 0.5772 | 0.2080 | 0.3693 | 0.5046 | 0.5091 | 0.0005 |
| 555 | 1.0137 | 0.9316 | 0.4276 | 31 | 55 | 0.6066 | 0.2056 | 0.4010 | 0.6044 | 0.2088 | 0.3956 | 0.5688 | 0.5525 | 0.0005 |
| 101 | 1.0173 | 0.9669 | 0.4407 | 34 | 58 | 0.6053 | 0.2047 | 0.4006 | 0.6014 | 0.2076 | 0.3938 | 0.5898 | 0.5739 | 0.0005 |
| 2718 | 1.0276 | 0.9882 | 0.4414 | 35 | 59 | 0.6135 | 0.2051 | 0.4084 | 0.6085 | 0.2111 | 0.3974 | 0.5464 | 0.5463 | 0.0005 |
| 77 | **1.0524** | 1.0129 | 0.4392 | 33 | 57 | 0.6267 | 0.2058 | 0.4210 | 0.6072 | 0.2072 | 0.4000 | 0.5739 | 0.5611 | 0.0005 |

### `d2_i::wsi_biology`, all twelve partitions

| split seed | **retention** | retention (global null) | site coverage | frequent sites | ind. S1 | ind. null | ind. excess | tr. S1 | tr. null | tr. excess | ind. S2 | tr. S2 | p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 31337 | **0.9221** | 0.8386 | 0.4870 | 39 | 0.5088 | 0.2058 | 0.3030 | 0.5370 | 0.2084 | 0.3285 | 0.4357 | 0.4578 | 0.0005 |
| **42** (published) | **0.9710** | 0.9022 | 0.4385 | 31 | 0.4898 | 0.2021 | 0.2877 | 0.5020 | 0.2056 | 0.2963 | 0.4170 | 0.4215 | 0.0005 |
| 77 | 0.9780 | 0.8928 | 0.4392 | 33 | 0.5119 | 0.2056 | 0.3063 | 0.5196 | 0.2064 | 0.3132 | 0.4258 | 0.4615 | 0.0005 |
| 23 | 0.9953 | 0.8834 | 0.4595 | 33 | 0.5041 | 0.2052 | 0.2989 | 0.5082 | 0.2079 | 0.3003 | 0.4071 | 0.4448 | 0.0005 |
| 101 | 1.0031 | 0.8990 | 0.4407 | 34 | 0.5059 | 0.2039 | 0.3020 | 0.5088 | 0.2077 | 0.3011 | 0.4610 | 0.4570 | 0.0005 |
| 11 | 1.0099 | 0.8995 | 0.4747 | 35 | 0.4902 | 0.2038 | 0.2864 | 0.4896 | 0.2060 | 0.2836 | 0.4328 | 0.4350 | 0.0005 |
| 4242 | 1.0230 | 0.9348 | 0.4711 | 38 | 0.5026 | 0.2051 | 0.2975 | 0.4981 | 0.2073 | 0.2908 | 0.3979 | 0.3886 | 0.0005 |
| 7 | 1.0255 | 0.9703 | 0.4356 | 32 | 0.4908 | 0.2034 | 0.2873 | 0.4862 | 0.2060 | 0.2801 | 0.3943 | 0.3937 | 0.0005 |
| 9 | 1.0406 | 0.9208 | 0.4457 | 33 | 0.4685 | 0.2069 | 0.2616 | 0.4582 | 0.2068 | 0.2514 | 0.3766 | 0.3541 | 0.0005 |
| 1000 | 1.0457 | 1.0014 | 0.4609 | 36 | 0.5485 | 0.2039 | 0.3447 | 0.5376 | 0.2079 | 0.3296 | 0.5433 | 0.5095 | 0.0005 |
| 2718 | 1.0516 | 0.9900 | 0.4414 | 35 | 0.4826 | 0.2061 | 0.2765 | 0.4716 | 0.2086 | 0.2630 | 0.4031 | 0.3644 | 0.0005 |
| 555 | **1.0794** | 1.0322 | 0.4276 | 31 | 0.5163 | 0.2057 | 0.3106 | 0.4953 | 0.2075 | 0.2877 | 0.4205 | 0.4352 | 0.0005 |

**Read the `d2_i` spread against the arms it is a ratio of.** `d2_i`'s channel is the smaller one —
excess 0.2616–0.3447 against `d2_h`'s 0.3640–0.4257 — so the same absolute partition-to-partition
wobble in the numerator and the denominator is a larger *fraction* of it. The inductive S1 itself
varies by only 0.080 across `d2_i`'s twelve partitions (0.4685–0.5485) and its matched control by
0.079 (0.4582–0.5376); the retention ratio inherits both. **This is a small-denominator effect, not
`d2_i` behaving differently out of sample**, and it is why the spread bar fires on the weaker artifact
first. It is reported as a limit on how precisely retention can be quoted at this cohort size, not as
evidence against the channel.

**The `none` arm, on every partition, for scale.** Retention of the unadjusted exposure block against
the same control: `d2_h` median **0.2602** (0.2211–0.3055), `d2_i` median **0.3521** (0.2991–0.3796).
Every inductive partition sits 3–4× above the best `none` partition. Whatever the spread is, it is
not the scale on which the adjustment question is decided.

## 4. Is the spread mechanistic? — the coverage question, answered two ways

### 4a. Within the twelve fixed-fraction partitions: no, and the correlation has the wrong sign

At f = 0.5 the discovery fold's site coverage moves little — **0.4276 to 0.4870** across twelve
partitions, a span of 5.9 percentage points — because coverage is mostly a property of TCGA's
site-size distribution and not of the draw. (It moves *more* than the ±3 pp §2 predicted: three of
the twelve sit 3.3–4.9 pp above the published partition's 43.85%, so that prediction is graded ✗
below. The direction is helpful — the seed sweep varies coverage somewhat more than I allowed for,
and retention still does not track it.) Spearman(retention, coverage)
= **−0.406 (p = 0.191)** on `d2_h` and **−0.259 (p = 0.417)** on `d2_i` — *negative*, i.e. the
opposite sign to the mechanism, and at n = 12 with essentially no power. `d2_i`'s two extremes make
the point without a correlation: the **worst** retention (0.9221) belongs to the partition with the
**best** coverage (0.4870, 39 frequent sites) and the **best** retention (1.0794) to the one with the
**worst** (0.4276, 31 sites). **Within the range a random partition can produce, the variation is not
coverage. It is sampling noise in a 1,382-row measurement, and no mechanism is invented for it.**

### 4b. When coverage is moved on purpose: yes, monotonically, and far too small to matter

`d2_h`, split seed 42, discovery fraction swept. Coverage is manipulated over 25 percentage points,
five times the range random partitions produce.

| discovery fraction | discovery rows | exposure rows | **site coverage** | frequent sites | operator design | **retention** | ind. S1 | ind. null | tr. S1 | tr. null | p | degeneracy flags |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **0.3** | 829 | 1,937 | **0.3015** | **15** | **39** | **0.9762** | 0.6089 | 0.1827 | 0.6128 | 0.1762 | 0.0005 | `few_frequent_sites`, `narrow_operator_design` |
| 0.5 | 1,384 | 1,382 | 0.4385 | 31 | 55 | 0.9966 | 0.6145 | 0.2067 | 0.6173 | 0.2080 | 0.0005 | — |
| **0.7** | 1,937 | 829 | **0.5476** | **53** | **77** | **1.0164** | 0.6164 | 0.2628 | 0.6158 | 0.2680 | 0.0005 | — |

**Monotone in coverage, and the slope is 0.0402 of retention across 24.6 points of coverage — about
0.0016 per point.** Prediction row 9 of the predeclaration ("f = 0.3 below f = 0.7 by more than
0.03") is **correct**: 1.0164 − 0.9762 = 0.0402. **So coverage does buy retention, and the effect is
real — and it is an order of magnitude too small to threaten the claim.** At f = 0.3 the operator's
design has collapsed from 55 columns to **39**, its frequent-site list from 31 sites to **15**, and
it can site-adjust only **30%** of the rows it scores — both §5.2 degeneracy flags fire, exactly as
they were written to — and the channel still reads retention **0.976** at the permutation floor. The
mechanism P4's 23:00 entry measured, where the same coverage loss costs the *site certificate* a
factor of 24, does not translate into the cross-block channel at any coverage this cohort can
produce.

## 5. The labels-only ceiling across partitions — the leg the verdict rests on

Share of the channel's own excess, `retention_of_excess(ceiling_adjusted, channel_arm)`, twelve
partitions, n = 1,382 in every cell.

| labels encoding | arm | `d2_h` median | `d2_h` range | sig. p < 0.05 | `d2_i` median | `d2_i` range | sig. |
|---|---|---:|---|---:|---:|---|---:|
| additive design, 108 cols | transductive | +0.2386 | **+0.0223 … +0.3777 (16.9×)** | 11/12 | +0.3091 | +0.0337 … +0.5065 (15.0×) | 11/12 |
| additive design | **inductive** | +0.2579 | +0.1792 … +0.4148 (2.3×) | 12/12 | +0.3722 | +0.2376 … +0.5492 (2.3×) | 12/12 |
| saturated cell design, 105 cols | transductive | +0.0998 | +0.0243 … +0.1920 | 8/12 | +0.1401 | +0.0368 … +0.2651 | 8/12 |
| saturated cell design | **inductive** | +0.1332 | +0.0213 … +0.4346 | 9/12 | +0.1758 | +0.0320 … +0.5297 | 9/12 |
| **operator's own frozen design** | transductive | +0.0456 | −0.0343 … +0.1678 | **3/12** | +0.0655 | −0.0449 … +0.2178 | **3/12** |
| **operator's own frozen design** | **inductive** | **−0.1221** | **−0.2343 … +0.0248** | **0/12** | **−0.1690** | −0.3221 … +0.0329 | **0/12** |

**Three things follow, and the third is the one that matters.**

1. **The 23:45 entry's central claim about the ceiling generalises: it does not rise out of sample.**
   The inductive arm's additive share is *lower* than its matched transductive control's on 4 of 12
   `d2_h` partitions (and 4 of 12 on `d2_i`) and higher on the other 8, medians 0.2579 vs 0.2386 —
   inside the partition-to-partition wobble of either, and nowhere near the rise the site-leakage
   mechanism predicted. The 23:45 entry's single-partition reading, where the inductive share was
   *below* the control on all three encodings, was one of the 4.
2. **The absolute ceiling share is the least stable quantity in this entry**, and it is *less* stable
   transductively (16.9×) than inductively (2.3×). Anyone quoting 6.0% or 11.2% is quoting one draw
   of a quantity that runs from 2.2% to 37.8% at fixed encoding and fixed n.
3. **On the encoding the adjustment spans, a confound-only representation carries nothing — on
   12 of 12 partitions, on both artifacts.** The frozen-design inductive ceiling is negative in 11 of
   12 (`d2_h`) and 11 of 12 (`d2_i`); the only positive readings are +0.0248 (`d2_h`) and +0.0329
   (`d2_i`), both on split seed 4242 and both at **p = 0.3688**; and **not one of the 24
   partition-arms is significant at 0.05.** The 23:45 entry called this
   "the strongest form of P1 §5's argument, and the first time it has been shown out of sample";
   §4/N5 of this predeclaration was written to kill it if it were partition-dependent. **It is not.**

## 6. Every §5 and §6 distrust check, discharged

| # | check | bar | measured over 24 partition-arms | fires? |
|---|---|---|---|:---:|
| §5.1 | the split actually changes | pairwise Jaccard < 0.9 | 0.3168–0.3450, 8/8 distinct, 12 distinct digests | no |
| §5.2 | frequent sites in the discovery fold | ≥ 20 | **31–39** on every f = 0.5 partition (15 at f = 0.3, flagged) | no |
| §5.2 | operator design width | ≥ 45 cols | **55–63** (39 at f = 0.3, flagged) | no |
| §5.2 | raw-vs-adjusted per-axis corr | 0.65–0.85 | `d2_h` 0.7531–0.7784, `d2_i` 0.7964–0.8183 | no |
| §5.2 | residual variance ratio | 0.50–0.70 | `d2_h` 0.5800–0.6200, `d2_i` 0.6686–0.6966 | no |
| §5.2 | axes with corr > 0.99 | 0 | **0** on every one of 24 arms | no |
| §5.3 | inductive null vs its control's | not < 0.70× | **0.9717–1.0008** — never more than 2.8% below | no |
| §5.4 | patients in both folds | 0 | **0**, asserted in the driver on all 25 runs | no |
| §6.1 | a partition that refuses | none | **none**; the replacement list {4242, 9, 77} was never needed for that purpose | no |
| §6.2 | is it the control that moved? | both arms tabulated | ind. S1 range 0.0600 vs control 0.0541 (`d2_h`), 0.0800 vs 0.0794 (`d2_i`) — the two arms move together | no |
| §6.3 | matched n | always | every comparison at n = 1,382 within its own partition | no |

**The adjustment is a real adjustment on every partition and is not the transductive one relabelled.**
`adjuster_agreement(inductive, transductive)` across twelve `d2_h` partitions: per-axis correlation
median **0.9548–0.9677**, relative Frobenius difference median **0.2862**,
`is_relabelled_incumbent` **False in 12/12** (and 12/12 on `d2_i`). S2 — directions fitted on one half
of the exposure fold and scored on the other, immune to S1's in-sample maximisation — moves by a
median of **−0.0022** (`d2_h`) and **−0.0008** (`d2_i`) between the arms, and is *higher* in the
inductive arm on 6 of 12 partitions on each artifact.

## 7. How the predictions did

| # | predicted (§2) | measured | |
|---|---|---|:---:|
| 1 | `d2_h` median retention 0.975 | **1.0102** | ✗ (too low by 0.035) |
| 2 | `d2_h` min 0.93 | **0.9872** | ✗ (too low) |
| 3 | `d2_h` max − min 0.07 | **0.0652** (0.0311 on the eight) | ✓ |
| 4 | all 8 `d2_h` in 0.90–1.05 (p = 0.60) | yes, 0.9966–1.0276 | ✓ |
| 5 | **not** all 8 `d2_h` in 0.97–1.00 (p = 0.75) | **correct — 6 of 8 sit above 1.00** | ✓ |
| 6 | `d2_i` median / min 0.960 / 0.90 | **1.0165 / 0.9221** | ✗ / ✓-ish |
| 7 | every arm at the 0.0005 floor (p = 0.90) | **24/24** | ✓ |
| 8 | coverage within ±3 pp of 43.8% (p = 0.80) | **42.76–48.70%**, i.e. −1.1 to **+4.9 pp**; 2 of the predeclared 8 (3 of 12) sit above +3 pp | ✗ |
| 9 | f = 0.3 below f = 0.7 by > 0.03 (p = 0.45) | **0.0402** | ✓ |
| 10 | frozen ceiling ≤ 0 or insignificant on all 8, inductive (p = 0.55) | **0/12 significant, both artifacts** | ✓ |

**Six of ten. Three of the four misses run the same way: I predicted the single-split figure was a
typical-to-lucky draw and that the distribution sat below it. It sits above it.** (The fourth,
row 8, is a miss in the direction that makes §4a's negative result stronger, not weaker.) Retention crosses
1.0 on 17 of 24 partitions, which I did not predict at all and which has a plain mechanism, measured
in §6: the inductive arm's within-cancer null median is systematically **0.3–2.8% *below*** its
matched in-sample control's (`null_ratio` 0.9717–1.0008, never above 1.001), while its S1 is
statistically indistinguishable — so the ratio of excesses lands just above one. §5.3 was written to
catch exactly this if it went the other way and by 30×; it went this way and by 2%.

## 8. What this changes

**The claim survives, its precision does not.** The morphology→molecular channel is intact under an
out-of-sample confound adjustment on **12 of 12** random partitions of the certification cohort, on
**both** artifacts, at the permutation floor in **24 of 24** arms, with the labels-only ceiling on the
adjustment's own design **insignificant in 24 of 24**. The narrowing this run was written to find —
"it survives *this* inductive adjustment" — **did not materialise**. What did materialise is that the
two quoted figures are one draw each from distributions of width 0.065 (`d2_h`) and 0.157 (`d2_i`),
and are near the bottom of both.

### Prose flagged, and deliberately NOT edited (multiple agents are in these files)

1. **`paper/P1_CALIBRA_DRAFT.md` §4.2, the residual-bound paragraph beginning "The channel itself,
   independent of the ceiling question, has now been validated out of sample."** (currently lines
   821–830; the sentence *"leaves it essentially unmoved: retention **0.9966** (`d2_h`) and **0.9710**
   (`d2_i`) against the transductive reading, on the identical 1,382 patients"*). This is the
   paragraph updated once today with the single-split result. **Both figures should become ranges over
   partitions**, e.g. *"…leaves it essentially unmoved on every one of twelve random discovery/exposure
   partitions of the certification cohort: retention 0.987–1.052 (`d2_h`, median 1.010) and 0.922–1.079
   (`d2_i`, median 1.017), with the single partition previously reported — 0.9966 and 0.9710 — the
   second-lowest of the twelve in each case, at p ≤ 0.0005 in every arm."* The `d2_i` range **must**
   be given as a range; a point quotation there is not supportable.
2. **Same file, §4.2 provenance footnote (currently lines 850–859), the sentence "retention reads
   0.997 under this paper's convention and 0.757 under a global-permutation convention; both are
   reported rather than one chosen silently."** The 0.757 is a **single-partition outlier**: across
   twelve partitions the global-null retention is 0.7572, then 0.9226 … 1.0129, median 0.9643, and the
   seed-42 value sits 0.165 below the next-worst. The sentence should say that the global-convention
   figure on the originally reported partition, 0.757, is the extreme minimum of twelve and that the
   median is 0.964 — otherwise the paper is offering its worst partition as the sceptic's reading.
3. **Same file, §4.2, the ceiling caveat sentence** (*"as a share of the channel's excess over that
   null they account for 11.2% and 6.0% respectively — when the labels are encoded in the same design
   columns the adjustment residualises against"*). The condition named there is correct and
   insufficient: **at fixed encoding and fixed n the share moves 16.9× across partitions** (additive,
   transductive, `d2_h`: 2.2% to 37.8%). The honest form is that the ceiling's absolute share is not a
   stable property of this cohort at all and is quotable only as "well below the channel" plus the
   frozen-design result, which *is* stable (0/24 significant).
4. **`NOTEBOOK_ENTRIES/inductive_channel_and_ceiling_result_20260804T2345Z.md` §12, first bullet
   ("One split. Discovery fraction 0.5 at seed 42 only")** — dischargeable, and should point here.
   Its §5's global-null discussion should carry the outlier finding from item 2.

**For P4.** Nothing here touches P4's condition 3. The site certificate's out-of-sample failure and
the channel's out-of-sample survival remain the two different claims about the same matrix that the
23:00 and 23:45 entries separated; this run adds only that the second one holds on every partition
and the first was never re-tested here.

## 9. Suite

Both readings are given because other agents' commits landed underneath this one while the runs were
in flight. `pytest morpheus/v2/tests morpheus/tests --ignore=morpheus/v2/tests/test_p2_figures.py -q`:

* at the commit the runs were launched from (`2b0ad53`, workspace `~/ws_p1sp`) →
  **623 passed, 0 failed in 69.04 s**;
* at this entry's final commit (workspace `~/ws_p1sp2`, freshly extracted) →
  **641 passed, 6 skipped, 0 failed in 72.08 s**. The extra 18 passes and 6 skips are the five
  commits other agents made in between (`6b3d8e7`…`e772316`, the P2 labelled-probe work); this run's
  code delta is unchanged between the two readings.

`test_p2_figures.py` run separately reads **1 passed, 27 errors in 2.36 s** at both commits, every
error `ModuleNotFoundError: No module named 'matplotlib'` — the known condition of `~/venv`.
**Nothing was installed into that environment.**

**This run's test delta is 12, not the full 13 by which the total exceeds the 23:45 entry's 610.**
`test_inductive_channel_split_stability.py` alone reads **12 passed**; the remaining +1 comes from
the four commits other agents made between `ce8f582` and `2b0ad53`
(`v2/tests/test_p2_{figures,floor_audit,labelled_probe,limit2_stress}.py` all moved), and is not
attributed here. `test_inductive_channel.py` (10), `test_inductive_adjustment.py` (20) and
`test_p4_inductive_wiring.py` (15) are unchanged and read **57 passed** together with the new file.

## 10. Files / provenance

Code change (one file, three default-preserving additions, none of which defines a statistic):
`v2/research/rebase/nature/p1_evidence/inductive_channel.py` — `--split-seed` threaded **only** into
`exposure_split` (`resolve_split_seed`, default `-1` → `--seed`, pinned by test); the P4
`prepare_state` bit-identity comparison recorded as *skipped with its reason* when the partition is
not P4's; and the exposure fold's site coverage read from `SitePooling.apply`'s existing report.
Tests `v2/tests/test_inductive_channel_split_stability.py`, 12 tests.
Probe `v2/research/rebase/nature/p1_evidence/split_overlap_probe.py` (imports `exposure_split` and
`_load_block`, defines nothing). Aggregation
`v2/research/rebase/nature/p1_evidence/split_stability_aggregate.py` — reads JSON and takes ranges;
the only computation in it that is not read straight out of a run is `scipy.stats.spearmanr` in §4a,
which is reported with its n and p and is not load-bearing for any conclusion.

Outputs `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/p1_evidence/split_stability/`
(26 run JSONs + `AGGREGATE.json` + `split_overlap_probe.json` + `logs/`), vendored to
`v2/research/rebase/nature/p1_evidence/out/split_stability/`.

Exact command, per run:

```
python -m morpheus.v2.research.rebase.nature.p1_evidence.inductive_channel \
  --artifact $P/runs/d2_final/artifacts/d2_{h,i}_seed42.npz \
  --targets $P/data/frozen_rna_targets.npz --output $OUT/<artifact>_<tag>.json \
  --state wsi_biology --partition test --n-jobs 5|6 \
  --labels-blocks additive_design saturated_cell_design frozen_discovery_design \
  [--split-seed S] [--discovery-fraction 0.3|0.7] [--skip-full]
```

## 11. Honest constraints on every number above

* **One cohort, one partition of it, `wsi_biology` only.** Twelve *sub*-partitions of the D2 `test`
  partition are twelve draws from **one** 2,766-patient cohort; they bound partition dependence, not
  cohort dependence. `full_biology` and `rna_biology` remain RNA-derived and near-circular at ~0.89
  and are still not a morphology→molecular measurement.
* **`min_site_count` left at the project default of 10 on every run.** The 23:45 entry named a
  sweep as "the obvious next experiment" and this run did not do it. The f = 0.3 / 0.7 arm moves the
  same underlying quantity — design width and site coverage — by a comparable amount, but it is not
  the same experiment.
* **The coverage gradient is one seed at three fractions**, not a fraction × seed grid. Its two
  non-published points are single runs and the monotonicity rests on three numbers.
* **The extension from 8 to 12 partitions was decided after seeing the eight.** It is reported as
  such, both figures are tabulated, and the verdict is graded on the predeclared eight.
* **Retention above 1.0 is not evidence that an out-of-sample adjustment is *better*.** It is the
  2%-lower inductive null of §7, and the difference in S1 between the arms is within the
  partition-to-partition wobble of either.
* **No number here is compared to a published n = 2,766 figure without its matched control between
  them**, and every retention is a ratio of excesses over each arm's own null, per
  `retention_of_excess`.
