"""TCGA-NSCLC vs ALCHEMIST: one instrument, pointed twice.

Runs, in this order and deliberately in this order:

1. the per-patient representation, built identically on both cohorts;
2. the frozen-artifact export in the ``v2/contracts.py`` NPZ shape;
3. **the cohort classifier** -- TCGA vs ALCHEMIST, plus the within-TCGA control, using
   ``v2/calibra/hest.py:cohort_classifier_auc`` unmodified, the same instrument that
   measured TCGA vs HEST at AUC 0.99999;
4. only then the channel, on each cohort separately, with its permutation null.

Step 3 comes before step 4 because a cohort-separability number that arrives after the
channel number is a number that can be argued with.

The representation is ``concat(mean, std)`` over a patient's raw 1536-d H-Optimus-0 tokens,
then PCA-256 within cohort.  Zero fitted parameters see any label.  This is the same
representation the dilution curve was measured on, which is the only reason the dilution
curve can be used to predict the penalty for tissue-level sampling.

The channel is never measured on the two cohorts pooled.  Each cohort is residualised and
correlated against its own targets, so a between-cohort direction -- however separable the
cohorts turn out to be in step 3 -- is not available to either measurement.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3].parent))

from morpheus.v2.baseline_exports import export_raw_hoptimus_baseline  # noqa: E402
from morpheus.v2.calibra.calibration import permutation_null  # noqa: E402
from morpheus.v2.calibra.hest import cohort_classifier_auc  # noqa: E402
from morpheus.v2.calibra.residualise import (  # noqa: E402
    confound_design, cross_fitted_residuals, pooled_tissue_source_site,
)
from morpheus.v2.calibra.spectral import top_canonical_correlation  # noqa: E402

TOKENS_PER_PATIENT = 30
POOL_MIN_COUNT = 10
# Measured null-corrected retained channel against foreign-tissue fraction, from
# v2/research/rebase/nature/DILUTION_LOWER_BOUND.md section 2.
DILUTION_D = (0.00, 0.10, 0.20, 0.30, 0.40, 0.60, 0.80)
DILUTION_RETAINED = (1.000, 0.999, 0.968, 0.905, 0.804, 0.607, 0.333)


# A DECLARED mapping from ALCHEMIST ICD-O histology onto the TCGA lung code space.
# Written down before the channel was measured, and reviewable as a list rather than a rule,
# because a cancer-label choice is exactly the class of decision that can set the answer
# without anyone noticing -- the coordinate-system work moved the whole D2 effect by
# +0.118-0.120 on a basis change alone.
#
# Adenosquamous carcinoma is genuinely both lineages, so it goes to OTHER rather than being
# forced to whichever of LUAD/LUSC would flatter the result.
ALCHEMIST_TO_TCGA_LUNG = {
    "Adenocarcinoma_NOS": "LUAD",
    "Adenocarcinoma_in_situ_NOS": "LUAD",
    "Acinar_adenocarcinoma": "LUAD",
    "Mucinous_adenocarcinoma": "LUAD",
    "Invasive_mucinous_adenocarcinoma": "LUAD",
    "Mixed_invasive_mucinous_and_non-mucinous_adenocarcinoma": "LUAD",
    "Lepidic_predominant_adenocarcinoma": "LUAD",
    "Adenocarcinoma_with_mixed_subtypes": "LUAD",
    "Squamous_cell_carcinoma_NOS": "LUSC",
}


def map_to_tcga_lung(histologies: np.ndarray) -> np.ndarray:
    return np.asarray([ALCHEMIST_TO_TCGA_LUNG.get(str(h), "OTHER") for h in histologies])


def pool_rare(labels: np.ndarray, minimum: int = POOL_MIN_COUNT) -> np.ndarray:
    """The `min_site_count` rule `pooled_tissue_source_site` uses, applied to any label."""
    counts = pd.Series(labels).value_counts()
    frequent = set(counts[counts >= minimum].index)
    return np.asarray([label if label in frequent else "OTHER" for label in labels])


def nominal_split(identifiers: np.ndarray) -> np.ndarray:
    """A deterministic 70/15/15 label, present only to satisfy the artifact contract.

    Nothing in this script fits on it: the representation has no fitted parameters and the
    target scoring is cohort-fit-free, so every measurement below uses all rows.  It is
    recorded as nominal rather than quietly reused as if it meant something.
    """
    out = []
    for identifier in identifiers:
        bucket = int(hashlib.sha256(str(identifier).encode()).hexdigest()[:8], 16) % 100
        out.append("train" if bucket < 70 else ("val" if bucket < 85 else "test"))
    return np.asarray(out)


def load_alchemist(staging: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[np.ndarray], dict]:
    ids, cancers, scanners, blocks = [], [], [], []
    skipped = []
    for path in sorted(glob.glob(str(Path(staging) / "*.npz"))):
        raw = np.load(path, allow_pickle=False)
        meta = json.loads(str(raw["meta"]))
        embeddings = np.asarray(raw["embeddings"], dtype=np.float32)
        if embeddings.shape[0] < TOKENS_PER_PATIENT or not np.isfinite(embeddings).all():
            skipped.append(meta["slide"])
            continue
        ids.append(meta["patient"])
        cancers.append(meta["cancer"])
        scanners.append(meta.get("slide_properties", {}).get("aperio.ScanScope ID", "NA"))
        blocks.append(embeddings[:TOKENS_PER_PATIENT])
    return (np.asarray(ids), np.asarray(cancers), np.asarray(scanners), blocks,
            {"skipped": skipped, "n_loaded": len(ids)})


def load_tcga(store_h5: str, metadata_parquet: str, patients: list[str],
              seed: int = 42) -> tuple[np.ndarray, list[np.ndarray]]:
    import h5py

    metadata = pd.read_parquet(metadata_parquet, columns=["patient_id", "slide_id", "row_idx"])
    metadata = metadata[metadata.patient_id.isin(set(patients))]
    rng = np.random.default_rng(seed)
    chosen_ids, chosen_rows = [], []
    for patient, group in metadata.groupby("patient_id", sort=True):
        rows = np.sort(group.row_idx.to_numpy())
        if len(rows) < TOKENS_PER_PATIENT:
            continue
        pick = np.sort(rng.choice(len(rows), size=TOKENS_PER_PATIENT, replace=False))
        chosen_ids.append(patient)
        chosen_rows.append(rows[pick])
    flat = np.concatenate(chosen_rows)
    order = np.argsort(flat)
    with h5py.File(store_h5, "r") as handle:
        values = handle["embeddings"][np.sort(flat)]
    restored = np.empty_like(values)
    restored[order] = values
    blocks = [restored[i * TOKENS_PER_PATIENT:(i + 1) * TOKENS_PER_PATIENT].astype(np.float32)
              for i in range(len(chosen_ids))]
    return np.asarray(chosen_ids), blocks


def pooled_representation(blocks: list[np.ndarray]) -> np.ndarray:
    return np.stack([np.concatenate([block.mean(axis=0), block.std(axis=0)]) for block in blocks])


def pca_reduce(matrix: np.ndarray, n_components: int, seed: int) -> np.ndarray:
    from sklearn.decomposition import PCA

    n = min(n_components, matrix.shape[0], matrix.shape[1])
    return PCA(n_components=n, random_state=seed).fit_transform(matrix).astype(np.float32)


def dilution_prediction(d: float) -> float:
    return float(np.interp(d, DILUTION_D, DILUTION_RETAINED))


def measure(name: str, representation: np.ndarray, targets: np.ndarray, cancers: np.ndarray,
            *, extra_columns: dict[str, np.ndarray] | None = None, n_components: int,
            n_permutations: int, seed: int, n_jobs: int) -> dict:
    frame = pd.DataFrame({"cancer": cancers})
    columns = ["cancer"]
    for key, values in (extra_columns or {}).items():
        frame[key] = values
        columns.append(key)
    design = confound_design(frame, columns)
    summary = permutation_null(representation, targets, design, strata=cancers,
                               n_permutations=n_permutations, n_components=n_components,
                               seed=seed, n_jobs=n_jobs)
    summary.update({
        "arm": name, "n_patients": int(representation.shape[0]),
        "n_targets": int(targets.shape[1]), "design_columns": columns,
        "n_design_columns": int(design.shape[1]),
        "n_cancer_levels": int(len(set(cancers.tolist()))),
        "unadjusted_top_cca": float(top_canonical_correlation(
            representation, targets, n_components=n_components)),
    })
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alchemist-staging", required=True)
    parser.add_argument("--alchemist-targets", required=True)
    parser.add_argument("--tcga-targets", required=True)
    parser.add_argument("--tcga-store-h5", required=True)
    parser.add_argument("--tcga-metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-components", type=int, default=16)
    parser.add_argument("--pca-components", type=int, default=256)
    parser.add_argument("--n-permutations", type=int, default=300)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--classifier-n", type=int, default=20000)
    parser.add_argument("--assumed-contamination", default="0.20,0.50")
    # In-sample top-CCA is capacity-inflated, and the inflation depends on n. ALCHEMIST has
    # 1,106 patients against TCGA-NSCLC's 841, so ALCHEMIST's permutation null sits LOWER and
    # its null-corrected excess is flattered relative to TCGA's by the size difference alone.
    # This subsamples ALCHEMIST to TCGA's n so the ratio is read at matched capacity.
    parser.add_argument("--alchemist-subsample", type=int, default=0,
                        help="subsample ALCHEMIST to this many patients (0 = use all)")
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    report: dict[str, object] = {"config": vars(args)}

    # ---- 1. representations ----------------------------------------------------------
    alch_targets = np.load(args.alchemist_targets, allow_pickle=False)
    tcga_targets = np.load(args.tcga_targets, allow_pickle=False)
    names = list(alch_targets["target_names"].astype(str))
    if names != list(tcga_targets["target_names"].astype(str)):
        raise ValueError("the two cohorts were not scored on the same target set")

    alch_ids, alch_cancers, alch_scanners, alch_blocks, alch_info = load_alchemist(args.alchemist_staging)
    report["alchemist_extraction"] = alch_info
    print(f"[alch] {len(alch_ids)} slides loaded, skipped={len(alch_info['skipped'])}", flush=True)

    tcga_patient_ids = list(tcga_targets["patient_ids"].astype(str))
    tcga_ids, tcga_blocks = load_tcga(args.tcga_store_h5, args.tcga_metadata, tcga_patient_ids,
                                      seed=args.seed)
    print(f"[tcga] {len(tcga_ids)} patients with >= {TOKENS_PER_PATIENT} tokens", flush=True)

    # Align each cohort to its own target table.  Cancer labels are taken FROM the target
    # block, not re-derived here, so the labels used by the confound design are provably the
    # ones the scores were filed under.
    def align(ids, blocks, extra, table):
        order = {identifier: index for index, identifier in enumerate(
            table["patient_ids"].astype(str))}
        keep = [index for index, identifier in enumerate(ids) if identifier in order]
        rows = np.asarray([order[ids[index]] for index in keep])
        return (np.asarray(ids)[keep], [blocks[index] for index in keep],
                np.asarray(table["scores"], dtype=np.float64)[rows],
                np.asarray(table["cancers"]).astype(str)[rows],
                {key: np.asarray(values)[keep] for key, values in extra.items()})

    alch_ids2, alch_blocks, alch_y, alch_cancers, alch_extra = align(
        alch_ids, alch_blocks, {"scanner": alch_scanners}, alch_targets)
    alch_cancers = pool_rare(alch_cancers)
    alch_scanners = pool_rare(alch_extra["scanner"])
    tcga_ids2, tcga_blocks, tcga_y, tcga_cancers, _ = align(
        tcga_ids, tcga_blocks, {}, tcga_targets)
    if len(set(tcga_cancers.tolist())) < 2:
        raise ValueError(f"TCGA arm collapsed to one cancer level: {set(tcga_cancers.tolist())}")
    print(f"[align] alchemist {len(alch_ids2)} x {alch_y.shape[1]} "
          f"({len(set(alch_cancers.tolist()))} histology levels), "
          f"tcga {len(tcga_ids2)} x {tcga_y.shape[1]} "
          f"({sorted(set(tcga_cancers.tolist()))})", flush=True)

    if args.alchemist_subsample and args.alchemist_subsample < len(alch_ids2):
        pick = np.sort(np.random.default_rng(args.seed).choice(
            len(alch_ids2), size=args.alchemist_subsample, replace=False))
        alch_ids2, alch_y = alch_ids2[pick], alch_y[pick]
        alch_cancers, alch_scanners = alch_cancers[pick], alch_scanners[pick]
        alch_blocks = [alch_blocks[i] for i in pick]
        print(f"[match-n] ALCHEMIST subsampled to {len(alch_ids2)} patients", flush=True)

    alch_raw = pooled_representation(alch_blocks)
    tcga_raw = pooled_representation(tcga_blocks)

    # ---- 2. artifacts in the consumable contract shape --------------------------------
    provenance = {"encoder": "bioptimus/H-optimus-0", "tokens_per_patient": TOKENS_PER_PATIENT,
                  "fov_microns": 128.0, "analysed_fov_microns_after_centre_crop": 112.0,
                  "output_px": 256, "jpeg_quality": 75, "jpeg_subsampling": "4:2:0"}
    export_raw_hoptimus_baseline(
        output / "alchemist_wsi_identity.npz", patient_ids=alch_ids2, cancers=alch_cancers,
        split=nominal_split(alch_ids2),
        pooled_mean=np.stack([b.mean(0) for b in alch_blocks]),
        pooled_std=np.stack([b.std(0) for b in alch_blocks]),
        source_provenance={**provenance, "cohort": "ALCHEMIST-ALCH",
                           "sampling": "whole-slide tissue mask, no tumour polygon exists",
                           "split_is_nominal": True})
    export_raw_hoptimus_baseline(
        output / "tcga_nsclc_wsi_identity.npz", patient_ids=tcga_ids2,
        cancers=tcga_cancers, split=nominal_split(tcga_ids2),
        pooled_mean=np.stack([b.mean(0) for b in tcga_blocks]),
        pooled_std=np.stack([b.std(0) for b in tcga_blocks]),
        source_provenance={**provenance, "cohort": "TCGA-LUAD+LUSC",
                           "sampling": "pathologist tumour polygons (TCGA-UT)",
                           "split_is_nominal": True})
    print("[artifacts] written", flush=True)

    # ---- 3. the mandatory control, BEFORE any channel number --------------------------
    rng = np.random.default_rng(args.seed)

    def flat_patches(blocks, n):
        stacked = np.concatenate(blocks, axis=0).astype(np.float64)
        if len(stacked) > n:
            stacked = stacked[rng.choice(len(stacked), size=n, replace=False)]
        return stacked / np.maximum(np.linalg.norm(stacked, axis=1, keepdims=True), 1e-8)

    n_each = min(args.classifier_n, len(alch_blocks) * TOKENS_PER_PATIENT,
                 len(tcga_blocks) * TOKENS_PER_PATIENT // 2)
    alch_patches = flat_patches(alch_blocks, n_each)
    tcga_all = np.concatenate(tcga_blocks, axis=0).astype(np.float64)
    pick = rng.choice(len(tcga_all), size=min(2 * n_each, len(tcga_all)), replace=False)
    tcga_pool = tcga_all[pick]
    tcga_pool = tcga_pool / np.maximum(np.linalg.norm(tcga_pool, axis=1, keepdims=True), 1e-8)
    cohort = {
        "n_per_cohort": int(n_each),
        "tcga_vs_alchemist_auc": float(cohort_classifier_auc(tcga_pool[:n_each], alch_patches, seed=0)),
        "tcga_vs_tcga_control_auc": float(cohort_classifier_auc(
            tcga_pool[:n_each], tcga_pool[n_each:2 * n_each], seed=0)),
        "normalisation": "row L2 on raw 1536-d H-Optimus-0 patch embeddings",
        "reference_tcga_vs_hest_auc": 0.99999,
    }
    report["cohort_classifier"] = cohort
    print(f"[COHORT CONTROL] TCGA vs ALCHEMIST AUC = {cohort['tcga_vs_alchemist_auc']:.5f}   "
          f"within-TCGA control = {cohort['tcga_vs_tcga_control_auc']:.5f}", flush=True)
    (output / "cohort_control.json").write_text(json.dumps(cohort, indent=2))

    # ---- 4. the channel --------------------------------------------------------------
    alch_x = pca_reduce(alch_raw, args.pca_components, args.seed)
    tcga_x = pca_reduce(tcga_raw, args.pca_components, args.seed)
    common = dict(n_components=args.n_components, n_permutations=args.n_permutations,
                  seed=args.seed, n_jobs=args.n_jobs)
    arms = {}
    arms["tcga_nsclc_cancer_only"] = measure(
        "tcga_nsclc_cancer_only", tcga_x, tcga_y, tcga_cancers, **common)
    print(f"[channel] tcga_nsclc observed={arms['tcga_nsclc_cancer_only']['observed_top_cca']:.4f} "
          f"null_median={arms['tcga_nsclc_cancer_only']['null_median']:.4f}", flush=True)
    arms["alchemist_cancer_only"] = measure(
        "alchemist_cancer_only", alch_x, alch_y, alch_cancers, **common)
    print(f"[channel] alchemist observed={arms['alchemist_cancer_only']['observed_top_cca']:.4f} "
          f"null_median={arms['alchemist_cancer_only']['null_median']:.4f}", flush=True)

    # Sensitivities: the site term TCGA has and ALCHEMIST does not, and the scanner id
    # which is the closest available analogue on the ALCHEMIST side.  Neither is the bar.
    tss, sites = pooled_tissue_source_site(tcga_ids2, min_site_count=POOL_MIN_COUNT)
    arms["tcga_nsclc_cancer_plus_site"] = measure(
        "tcga_nsclc_cancer_plus_site", tcga_x, tcga_y, tcga_cancers,
        extra_columns={"tss": tss}, **common)
    arms["alchemist_cancer_plus_scanner"] = measure(
        "alchemist_cancer_plus_scanner", alch_x, alch_y, alch_cancers,
        extra_columns={"scanner": alch_scanners}, **common)

    # Label-granularity sensitivity. The primary ALCHEMIST arm adjusts on ALCHEMIST's own
    # pooled ICD-O histology; this one collapses it onto the TCGA lung code space
    # {LUAD, LUSC, OTHER} via the declared mapping above. If the channel moves between the
    # two, the cancer-label choice is setting the answer and must be reported as such.
    alch_mapped = map_to_tcga_lung(alch_cancers)
    arms["alchemist_mapped_to_tcga_lung_codes"] = measure(
        "alchemist_mapped_to_tcga_lung_codes", alch_x, alch_y, alch_mapped, **common)
    report["label_mapping"] = {
        "mapping": ALCHEMIST_TO_TCGA_LUNG,
        "levels_after_mapping": {level: int((alch_mapped == level).sum())
                                 for level in sorted(set(alch_mapped.tolist()))},
        "primary_arm_uses": "ALCHEMIST's own pooled ICD-O histology, fitted within ALCHEMIST",
        "why_not_the_persisted_tcga_operator":
            "runs_misc/tcga_operators/*cancer_only* was fitted on the 2,530-row TCGA TEST "
            "split, whose 22 design columns are its 21 cancers plus cancer_nan and contain "
            "NEITHER LUAD NOR LUSC (the test split holds 0 of each). TCGA's own NSCLC "
            "patients raise UnseenLevelError against it exactly as ALCHEMIST does, so it "
            "cannot adjust either arm of a lung-versus-lung comparison.",
    }

    # ---- 5. the predeclared bar ------------------------------------------------------
    reference = arms["tcga_nsclc_cancer_only"]
    candidate = arms["alchemist_cancer_only"]
    denominator = reference["excess_over_null_median"]
    ratio = float(candidate["excess_over_null_median"] / denominator) if denominator > 0 else float("nan")
    low, high = (float(v) for v in args.assumed_contamination.split(","))
    predicted = (dilution_prediction(high), dilution_prediction(low))
    p = candidate["permutation_p"]
    if p < 0.01 and ratio >= 0.60:
        verdict = "REPLICATES"
    elif p < 0.01 and ratio >= 0.30:
        verdict = "ATTENUATED BUT PRESENT"
    else:
        verdict = "FAILS TO REPLICATE"
    report["arms"] = arms
    report["verdict"] = {
        "R": ratio, "permutation_p": p, "verdict": verdict,
        "assumed_contamination_band": [low, high],
        "dilution_predicted_retained_band": list(predicted),
        "bar": {"REPLICATES": "p<0.01 and R>=0.60",
                "ATTENUATED BUT PRESENT": "p<0.01 and 0.30<=R<0.60",
                "FAILS TO REPLICATE": "p>=0.01 or R<0.30"},
        "predeclared_in": "NOTEBOOK_ENTRIES/PREDECLARED_alchemist_external_replication_20260804T1830Z.md",
    }
    mapped = arms["alchemist_mapped_to_tcga_lung_codes"]
    report["verdict"]["R_under_mapped_labels"] = (
        float(mapped["excess_over_null_median"] / denominator) if denominator > 0 else float("nan"))
    # Written last, so everything computed above is actually in the artifact. The first run
    # dumped the JSON before this line and the label-sensitivity ratio existed only in stdout.
    (output / "alchemist_channel.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\n[VERDICT] R={ratio:.3f}  p={p:.4f}  ->  {verdict}", flush=True)
    print(f"[label sensitivity] R under the declared LUAD/LUSC/OTHER mapping = "
          f"{report['verdict']['R_under_mapped_labels']:.3f} "
          f"(primary, ALCHEMIST's own histology: {ratio:.3f})", flush=True)
    print(f"[VERDICT] dilution predicted retained band for d in [{low},{high}] = "
          f"{predicted[0]:.3f}..{predicted[1]:.3f}", flush=True)


if __name__ == "__main__":
    main()
