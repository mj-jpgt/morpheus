## P1 §4.2's twelve numbers come from `runs/d2_final/artifacts/`, SHA-256 `4a18b94f…` and `028e8635…` — established three independent ways; the cancer-type pair `0.463 → 0.035` comes from an artifact we cannot identify, and that is now written into the papers

**Logged:** 2026-08-04 22:10 UTC.
**Question:** three files are called `d2_h_seed42.npz`, with three hashes and three raw joint-LDA
values (0.3633 / 0.1782 / 0.3785). The paper names neither path nor hash. Which artifact did each
published number come from?
**How obtained:** box `150.136.45.194`, workspace `~/ws_ind` from
`git -c core.autocrlf=false archive HEAD`, per-file blob SHA-1 verified. `~/venv`, threads capped to 1.
Every statistic imported from `v2/calibra/confound_certificate.py` (`lda_oof_balanced_accuracy`,
`_stratified_folds`, `_encode`) and `v2/calibra/residualise.py` (`pooled_tissue_source_site`); nothing
computed inline.

---

### 1. The three copies, hashed

`find` over `/home/ubuntu`, `/lambda` and `/mnt` returns exactly three. (`~/e0_run` is a symlink into
`/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/e0_run`, so the `~/…` and `/lambda/…`
forms are the same bytes, not extra copies.) There are **zero** copies in the checkout.

| path (under `/lambda/nfs/geeg/biorag3_persistent_20260711/morpheus_phase_d/`) | size | mtime (UTC) | SHA-256 |
|---|---:|---|---|
| `runs/d2_final/artifacts/d2_h_seed42.npz` | 18,348,420 | 2026-08-01 20:33:24 | `4a18b94f1017b85dd576f30ee8e3caf92d7897630a7054efb70166191cbe69e3` |
| `e0_run/d2_v3/recovered_artifacts/d2_h_seed42.npz` | 18,343,374 | 2026-08-03 01:11:46 | `053490d685bf0dc47f2094831048db2bb884fe99f7ada3b57508ca23b561b899` |
| `e0_run/d2_v3/d2_v3_s42/artifacts/d2_h_seed42.npz` | 18,333,599 | 2026-08-03 07:57:39 | `e81f4496f82c503a0dd1833e77cde2ea383cf79b0a6a7423a11c977c7f8f2625` |

`d2_i_seed42.npz` differs too: `028e8635465dd3c6…` (`d2_final`) against `b49dc3efaf0a25dd…`
(`d2_v3_s42`); `recovered_artifacts/` has no plain `d2_i_seed42.npz` at all, only
`d2_i_seed42_EPOCH13.npz`.

### 2. Which artifact each published number came from

**All twelve cells of the §4.2 table come from `runs/d2_final/artifacts/`**, established three ways
that are independent of each other:

1. **The run script names it.** `p1_evidence/run_track1.sh:6` sets
   `D2=…/runs/d2_final/artifacts` and passes `$D2/d2_h_seed42.npz $D2/d2_i_seed42.npz` to
   `python -m morpheus.v2.calibra.confound_certificate`, with and without `--residualise`. This is a
   record, not an inference.
2. **The run log prints the table.** `p1_evidence/logs/t1_cert_raw.log` contains
   `[d2_h_seed42::wsi_biology] axes=256 breach=17 per_axis_max=0.0532 chance=0.0118 joint=0.3633`,
   i.e. §4.2's row verbatim, and five more matching the other five states.
3. **Re-running the estimator discriminates.** Raw joint LDA recomputed over all three copies, all six
   states, `--partition test`, seed 42, 5 folds:

| copy | d2_h wsi | d2_h full | d2_h rna | d2_i wsi | d2_i full | d2_i rna |
|---|---:|---:|---:|---:|---:|---:|
| **`runs/d2_final/artifacts/`** | **0.363263** | **0.263038** | **0.256303** | **0.234848** | **0.268867** | **0.274354** |
| published in §4.2 | 0.3633 ✓ | 0.2630 ✓ | 0.2563 ✓ | 0.2348 ✓ | 0.2689 ✓ | 0.2744 ✓ |
| `e0_run/d2_v3/d2_v3_s42/artifacts/` | 0.378527 | 0.247556 | 0.253456 | 0.317895 | 0.277138 | 0.274867 |
| `e0_run/d2_v3/recovered_artifacts/` | 0.178154 | 0.253576 | 0.254231 | — | — | — |

**Six of six match `d2_final` to four decimals and none of the others do.** This also independently
reproduces the 0.1782 and 0.3785 figures, which until now existed only in the 21:00 entry's §6 table
and in no file on disk.

### 3. Why they differ

The `d2_v3` set is the **August re-run** (2026-08-03) and `d2_final` (2026-08-01) predates it. From
the `.diagnostics.json` sidecars, the three were written from **three different commits** —
`503c36b6` (dirty), `fa978460` (dirty), `3853d40d` (clean) — sharing one `configuration_sha256`
(`767fb363…`). They are different artifacts that share a filename, not copies of one artifact, which
is what the project's own recorded finding that training on this stack is not seed-reproducible
predicts. The gulf is structural, not cosmetic: effective rank 20.14 (`d2_v3_s42`) against 7.92
(`recovered_artifacts`).

Two further pointer defects found on the way, recorded because they were found:

* `runs/d2_final/artifacts/d2_h_seed42.npz.diagnostics.json` records its own `artifact` field as
  `/home/ubuntu/e0_run/d2_final/artifacts/d2_h_seed42.npz` — **a path that no longer exists**.
* §4.2 cited its command as `confound_certificate …` as though a console script existed. There is no
  `pyproject.toml` or `setup.py` in the repository, so there is none; the real invocation is
  `python -m morpheus.v2.calibra.confound_certificate`. Corrected.

### 4. One published number cannot be identified, and the papers now say so

**`0.463 → 0.035` (cancer-type balanced accuracy, chance 0.048, n = 2,530) has no identified
artifact.** `v2/research/rebase/nature/PHASE1_RESULT.md` states it in prose and names no path and no
hash. No run output under `p1_evidence/`, `p1_out/` or `e0_run/` records a cancer-type balanced
accuracy of 0.463 or 0.035; an exhaustive grep returns only unrelated coincidences (a dilution effect
size `0.1782483…`, a k-NN accuracy `0.1782898…`, a training loss `-0.3785374…`). The cohort is a
different one from the site arm — n = 2,530 against n = 2,766 — so it cannot be inherited from it.

It is **not** attributed to whichever artifact reproduces it. It is marked `**artifact not
identified**` at all four places P1 quotes it, carries a `PROVENANCE UNRESOLVED` block in §4.2 and at
its origin in `PHASE1_RESULT.md`, is flagged on P1_FIGURES panel (d), and is pinned with
`"artifact": None` plus a written reason in the new test. It must be regenerated against a hashed
artifact or withdrawn before submission.

### 5. What changed

**Papers now name path and SHA-256 for every quoted number.**
`paper/P1_CALIBRA_DRAFT.md` §4.2 gains an *Artifacts, by content hash* table for the two published
artifacts plus the three-copy discrimination table, and the corrected command line; the four other
places quoting 0.3633 (abstract ×2, §1.4, §6) carry an inline
``(artifact `runs/d2_final/artifacts/d2_h_seed42.npz`, SHA-256 `4a18b94f1017b85d…`; §4.2)``.
`paper/P1_FIGURES.md` F2 gains the same hash table before its panel list, and panel (d) is marked.
`v2/research/rebase/nature/TRACK1_NEGATIVE_CONTROLS.md` §T1.3 gains the four-row hash table.
`v2/research/rebase/nature/PHASE1_RESULT.md` gains the `PROVENANCE UNRESOLVED` block at the origin.

**New test — `v2/tests/test_paper_artifact_digests.py`** (10 tests, all passing). The existing
`test_paper_paths_resolve.py` checks that a cited path *exists*, which is why it could never have
caught this: §4.2 cited only its output directory, `p1_evidence` is in `BOX_TREES` and so skipped
outright, and the artifact was never named at all. The sibling checks identity instead of existence:

* `test_every_artifact_basename_in_a_draft_is_accompanied_by_a_hash` — every mention of
  `d2_h_seed42.npz` / `d2_i_seed42.npz` must have a SHA-256 (full or 16-hex) within 1,600 characters.
  A bare filename is not an identifier when three files share it.
* `test_every_pinned_number_is_hash_identified_in_its_paper` — the full digest must appear in the
  document *and* a digest must sit within 1,600 characters of every occurrence of the number.
* `test_an_unidentified_number_is_declared_unidentified` — a number with `"artifact": None` must be
  marked as such near every quotation. This is the check that keeps the registry honest: the failure
  mode worth guarding is not an absent hash, it is a *confident wrong* one.
* `test_recorded_digests_match_the_box_when_reachable` — re-hashes the real files when the NFS mount
  is present, skips with an explanatory message in a checkout.
* plus registry well-formedness and anti-vacuity checks.

All five hashes above are pinned in `ARTIFACTS`, including the two that are **not** published,
so a future edit citing them is caught rather than accepted. `test_paper_paths_resolve.py` gains the
two basenames to `BOX_OUTPUT_BASENAMES` with a reason pointing at the stricter rule.

### Honest constraints

The three-way identification is for the site certificate only. The §4.2 *adjusted* column was
confirmed against the persisted `certificate_adjusted/confound_certificate.json` rather than re-run
from the artifacts, because reproducing it needs the 1,000-permutation run; the raw column, which is
deterministic, was re-run. The scope of the new test is `PINNED_NUMBERS` — three numbers today, not
every number in P1; it is a mechanism plus a seed registry, and every further number added to the
registry is one more that cannot be quoted without its bytes. The unresolved `0.463 → 0.035` is
recorded, not fixed: regenerating it needs the Phase 1 cancer-type check re-run against a hashed
artifact, which was not in scope here.
