# P2 labelled linear probe — run output, vendored 2026-08-05

The measurement `paper/P2_RANK_DRAFT.md` §2.5 and §6.2 both marked as never run: the reference
standard RankMe and LiDAR were validated against, fitted on the same twelve frozen artifacts, the
same three co-trained views and the same cohort rule §4.1–§4.7 use, plus the five same-seed retrains
§4.1's floors are measured on.

* **Predeclaration:** `NOTEBOOK_ENTRIES/PREDECLARED_p2_labelled_linear_probe_20260805T0040Z.md`
  (`6b3d8e7`), committed before any statistic here existed.
* **Result:** `NOTEBOOK_ENTRIES/p2_labelled_linear_probe_result_20260805T0150Z.md`.
* **Runner:** `../p2_labelled_probe.py`; tests `v2/tests/test_p2_labelled_probe.py`.

These files are here because the paper's own standard is that every quoted number traces to a file in
the repository. `MANIFEST.json` carries a SHA-256 for each.

| file | what it is |
|---|---|
| `P2_LABELLED_PROBE.json` | the merged run + `_summary` (floors, per-pair table, agreement counts, channel/probe conflicts). **This is the file the notebook entry quotes.** |
| `shard_s1..s4.json` | the four shards it was merged from — s1 = D2 s42/s43, s2 = D2 s44 + D1 s42, s3 = D1 s43/s44, s4 = the five `d1_envelope` retrains |
| `shard_reps_panel.json` | a **redundant** second pass over the five retrains. The panel endpoints were already scored by `s4`; `merge()`'s duplicate guard refused it, which is how that was discovered. Kept because it reproduces `s4`'s panel values bit-for-bit and is therefore a determinism check on the estimator. |

## How it was produced

Box `150.136.45.194`, `~/ws_p2probe`, workspace verified **753/753 tracked files byte-equal to HEAD**
by per-file git blob SHA-1 before anything ran. `~/venv`, threads capped to 1 per
`NOTEBOOK_ENTRIES/operational_shared_box_rules_20260804T0730Z.md`. CPU only; no GPU and no retraining
— every probe is fitted on frozen exported embeddings.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYTHONPATH=$WS python -m morpheus.v2.research.rebase.p2.p2_labelled_probe \
  --labels .../e1_endpoints/inputs/e1_endpoint_labels.parquet \
  --targets ~/e0_run/data/frozen_rna_targets.npz \
  --views wsi_biology rna_biology full_biology --null-views wsi_biology \
  --endpoints mut_TP53 \
  --wsi-only-endpoints grade_high stage_late mut_ATM mut_KMT2D mut_ARID1A \
  --n-permutations-logistic 20 --n-permutations-lda 200 --n-boot 1000 --n-permutations-auroc 1000 \
  --output out/shard_sN.json --artifacts LABEL=path ...

python -m morpheus.v2.research.rebase.p2.p2_labelled_probe \
  --merge out/shard_s1.json out/shard_s2.json out/shard_s3.json out/shard_s4.json \
  --output out/P2_LABELLED_PROBE.json
```

Sharding is safe only because artifacts are independent until `summarise()` runs; the floors and the
agreement table are **not** shardable, since a floor built from a subset of the five retrains is not
the floor. `--merge` exists so that sharding goes through the same `summarise()` a single-process run
uses, and the test pins that a three-way shard reproduces the single-process summary exactly.

## Reproduction check

Every published number this run recomputes reproduces to every published digit: all twelve artifacts'
canonical R1 and untrained-40 channel (§4.1, §4.6), the five retrains' rank and channel, §4.1a's three
per-view floors (3.295× / 1.019× / 1.020×) and the 1.055× channel spread, and E1's six `d2_h_seed42`
endpoint AUROCs.
