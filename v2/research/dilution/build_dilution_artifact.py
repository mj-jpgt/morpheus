"""Build tissue-dilution CALIBRA artifacts from tumour + normal H-Optimus tokens.

Reproduces the ``raw_hoptimus_meanstd`` baseline representation -- the unweighted
``concat(mean, std)`` global skip path of ``v2/model.py:130-165`` -- over bags
that have been contaminated with a controlled fraction of foreign tissue.

Dilution level ``d`` is the fraction of the FINAL bag that is not the patient's
own tumour, so a patient with ``T`` tumour patches receives
``round(d/(1-d) * T)`` foreign patches.  Levels are NESTED: the patches used at
10% are a prefix of those used at 60%, so the curve is not confounded by
independent draw noise between levels.

Foreign patches are drawn donor-slide-first, not patch-first: a patient takes as
many patches as it needs from one randomly assigned donor slide before moving to
the next.  That mimics a whole-slide bag which includes a contiguous stretch of
the wrong tissue, rather than an unrealistic mosaic of one patch from each of
forty donors.

FOUR ARMS, because the obvious experiment has a preparation confound.  The
tumour store is 100% FFPE diagnostic (``DX``) slides; 98.5% of TCGA's public
solid-tissue-normal WSIs are frozen sections (``TS``/``BS``).  Frozen and FFPE
tissue are trivially separable in H-Optimus space, so a normal-tissue arm alone
cannot distinguish "non-tumour biology" from "different slide preparation".

  ``pooled``          pan-cancer frozen normal.  Upper bound: biology + preparation.
  ``matched``         organ-matched frozen normal, on the test cancers with
                      enough normal donors.
  ``foreign_tumour``  same-cancer, different-patient tumour patches drawn from
                      the SAME store.  Preparation-matched and pipeline-matched:
                      isolates the cost of adding patches that carry no
                      information about this patient from the domain shift.
  ``dx_normal``       the 14 FFPE diagnostic normal slides that do exist.
                      Preparation-matched non-tumour, but a very thin donor pool.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ARTIFACT_VERSION = 2


def _normalise(value: np.ndarray) -> np.ndarray:
    """Row L2 normalisation, identical to ``v2/baseline_exports._normalise``."""
    value = np.asarray(value, dtype=np.float32)
    return value / np.maximum(np.linalg.norm(value, axis=1, keepdims=True), 1e-8)


def _write(output: Path, *, patient_ids, cancers, split, states: dict[str, np.ndarray],
           method: str, config: dict) -> Path:
    trained = np.asarray(tuple(dict.fromkeys(sorted(states))), dtype=str)
    manifest = {"artifact_version": ARTIFACT_VERSION, "method": method,
                "trained_states": trained.tolist(),
                "fit_population": "none_no_fitted_representation", "config": config,
                "cohort_digest": sha256(json.dumps(
                    {"patient_ids": list(map(str, patient_ids)), "cancers": list(map(str, cancers)),
                     "split": list(map(str, split))}, sort_keys=True).encode()).hexdigest()}
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, patient_ids=np.asarray(patient_ids, dtype=str),
                        cancers=np.asarray(cancers, dtype=str), split=np.asarray(split, dtype=str),
                        trained_states=trained, artifact_version=np.asarray(ARTIFACT_VERSION),
                        manifest_json=np.asarray(json.dumps(manifest, sort_keys=True)), **states)
    return output


def load_normal_pool(staging: Path) -> tuple[np.ndarray, pd.DataFrame]:
    """Concatenate per-slide normal embeddings into one pool plus a donor index."""
    blocks, rows = [], []
    offset = 0
    for path in sorted(staging.glob("*.npz")):
        data = np.load(path, allow_pickle=True)
        meta = json.loads(str(data["meta"]))
        vectors = np.asarray(data["embeddings"], dtype=np.float32)
        if vectors.ndim != 2 or vectors.shape[1] != 1536 or not np.isfinite(vectors).all():
            raise ValueError(f"invalid normal embedding block in {path}")
        blocks.append(vectors)
        rows.append({"slide": meta["slide"], "patient": meta["patient"], "cancer": meta["cancer"],
                     "preparation": meta["slide"].split("-")[5][:2], "mpp": meta["mpp"],
                     "crop_px": meta["crop_px"], "start": offset, "n": len(vectors)})
        offset += len(vectors)
    if not blocks:
        raise RuntimeError(f"no normal embeddings under {staging}")
    return np.concatenate(blocks, axis=0), pd.DataFrame(rows)


def load_tumour_tokens(embeddings_h5: Path, metadata_parquet: Path, patient_ids, cancers):
    """Return the cohort's tumour tokens, per-patient spans and per-slide spans."""
    metadata = pd.read_parquet(metadata_parquet, columns=["patient_id", "slide_id", "row_idx"])
    wanted = {str(p): str(c) for p, c in zip(patient_ids, cancers)}
    metadata = metadata[metadata.patient_id.astype(str).isin(wanted)].sort_values(["patient_id", "slide_id", "row_idx"])
    with h5py.File(embeddings_h5, "r") as handle:
        tokens = np.asarray(handle["embeddings"][:], dtype=np.float32)[metadata.row_idx.to_numpy()]
    metadata = metadata.reset_index(drop=True)
    metadata["local"] = np.arange(len(metadata))
    bags, slide_rows = {}, []
    for patient, group in metadata.groupby("patient_id", sort=False):
        bags[str(patient)] = tokens[group.local.to_numpy()]
        for slide, sub in group.groupby("slide_id", sort=True):
            local = sub.local.to_numpy()
            slide_rows.append({"slide": str(slide), "patient": str(patient),
                               "cancer": wanted[str(patient)], "preparation": "DX",
                               "start": int(local[0]), "n": int(len(local))})
    missing = set(wanted) - set(bags)
    if missing:
        raise RuntimeError(f"{len(missing)} cohort patients have no tumour tokens, e.g. {sorted(missing)[:3]}")
    frame = pd.DataFrame(slide_rows)
    if not (np.diff(np.sort(np.concatenate([np.arange(r.start, r.start + r.n) for r in frame.itertuples()]))) == 1).all():
        raise RuntimeError("tumour slide spans are not a contiguous partition of the token block")
    return tokens, bags, frame


def patient_rng(seed: int, patient: str) -> np.random.Generator:
    """Per-patient RNG keyed by the patient ID, not by row position.

    Keyed by identity so the SAME patient draws the SAME foreign patches in the
    mean/std arm and in the trained-model arm, even though those two cohorts are
    ordered differently and are not the same size.
    """
    digest = sha256(patient.encode()).digest()[:8]
    return np.random.default_rng([seed, int.from_bytes(digest, "big")])


def normal_draw_order(donor_rows: pd.DataFrame, rng: np.random.Generator, needed: int) -> np.ndarray:
    """Row indices for one patient's foreign patches, donor-slide-first."""
    if needed <= 0 or donor_rows.empty:
        return np.zeros(0, dtype=np.int64)
    order = rng.permutation(len(donor_rows))
    picked: list[np.ndarray] = []
    total = 0
    for position in order:
        row = donor_rows.iloc[int(position)]
        picked.append(int(row.start) + rng.permutation(int(row.n)))
        total += int(row.n)
        if total >= needed:
            break
    return np.concatenate(picked)[:needed]


def build_arm(*, name: str, patient_ids, cancers, split, bags, pool, donors, levels,
              donor_lookup, seed: int, output: Path, config_extra: dict) -> dict:
    states: dict[str, np.ndarray] = {}
    per_level: dict[str, dict] = {}
    test = np.asarray(split) == "test"
    for level in levels:
        mean = np.zeros((len(patient_ids), 1536), dtype=np.float32)
        std = np.zeros((len(patient_ids), 1536), dtype=np.float32)
        drawn = np.zeros(len(patient_ids), dtype=np.int64)
        tumour_n = np.zeros(len(patient_ids), dtype=np.int64)
        for index, patient in enumerate(patient_ids):
            tumour = bags[str(patient)]
            tumour_n[index] = len(tumour)
            needed = int(round(level / (1.0 - level) * len(tumour))) if level > 0 else 0
            rng = patient_rng(seed, str(patient))
            rows = normal_draw_order(donor_lookup(str(patient), str(cancers[index])), rng, needed)
            drawn[index] = len(rows)
            bag = tumour if len(rows) == 0 else np.concatenate([tumour, pool[rows]], axis=0)
            mean[index] = bag.mean(axis=0, dtype=np.float64)
            std[index] = bag.std(axis=0, dtype=np.float64)
        key = f"wsi_dilution_{int(round(level * 100)):03d}"
        states[key] = _normalise(np.concatenate([mean, std], axis=1))
        achieved = drawn[test] / np.maximum(drawn[test] + tumour_n[test], 1)
        per_level[key] = {"level": level,
                          "achieved_foreign_fraction_test_median": float(np.median(achieved)),
                          "achieved_foreign_fraction_test_min": float(achieved.min()) if test.any() else 0.0,
                          "foreign_patches_total": int(drawn.sum()),
                          "foreign_patches_test": int(drawn[test].sum())}
        print(f"[{name}] {key}: foreign patches (test) = {per_level[key]['foreign_patches_test']}, "
              f"achieved median fraction = {per_level[key]['achieved_foreign_fraction_test_median']:.3f}", flush=True)
    config = {"arm": name, "pooling": "mean_plus_std_unweighted_global_skip_path",
              "dilution_levels": list(levels),
              "dilution_definition": "fraction_of_final_bag_that_is_not_this_patients_tumour",
              "nested_levels": True, "donor_assignment": "slide_first",
              "n_donor_slides": int(len(donors)), "n_donor_patches": int(len(pool)),
              "donor_preparation": donors.preparation.value_counts().to_dict(),
              "seed": seed, "per_level": per_level, **config_extra}
    _write(output, patient_ids=patient_ids, cancers=cancers, split=split, states=states,
           method=f"dilution_{name}", config=config)
    return config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-artifact", required=True,
                        help="npz supplying the exact patient_ids/cancers/split of the reference cohort")
    parser.add_argument("--embeddings-h5", required=True)
    parser.add_argument("--metadata-parquet", required=True)
    # Not required: the `foreign_tumour` arm needs no normal tissue at all.
    parser.add_argument("--normal-staging", default=None)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--levels", default="0.0,0.10,0.20,0.40,0.60")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-matched-donors", type=int, default=10)
    parser.add_argument("--arms", default="pooled,matched,foreign_tumour,dx_normal")
    args = parser.parse_args()

    levels = tuple(float(v) for v in args.levels.split(","))
    arms = [a for a in args.arms.split(",") if a]
    cohort = np.load(args.cohort_artifact, allow_pickle=True)
    patient_ids = np.asarray([str(p) for p in cohort["patient_ids"]])
    cancers = np.asarray([str(c) for c in cohort["cancers"]])
    split = np.asarray([str(s) for s in cohort["split"]])

    # `foreign_tumour` draws entirely from the tumour store, so it must not
    # require a normal-tissue staging directory. Loading the pool unconditionally
    # made the one GPU-free, no-new-data arm unrunnable on its own.
    needs_normals = bool({"pooled", "matched", "dx_normal"} & set(arms))
    if needs_normals:
        all_pool, all_donors = load_normal_pool(Path(args.normal_staging))
        # One staging directory holds both preparations; the arms are separated by
        # the slide-type field of the barcode, never by which directory it came from.
        frozen = all_donors.preparation.isin(("TS", "BS"))
        pool, donors = all_pool, all_donors[frozen]
        dx_donors = all_donors[all_donors.preparation == "DX"]
    else:
        all_pool = pool = np.empty((0, 0), dtype=np.float32)
        all_donors = donors = dx_donors = pd.DataFrame(
            columns=["patient", "preparation", "cancer", "slide", "start", "stop"])
    tumour_tokens, bags, tumour_slides = load_tumour_tokens(
        Path(args.embeddings_h5), Path(args.metadata_parquet), patient_ids, cancers)
    print(f"[normal pool] {len(all_donors)} slides ({all_donors.preparation.value_counts().to_dict()}), "
          f"{len(all_pool)} patches, {all_donors.patient.nunique()} donors; "
          f"frozen {len(donors)} / DX {len(dx_donors)}", flush=True)
    print(f"[cohort] {len(patient_ids)} patients, {(split == 'test').sum()} test, "
          f"{len(tumour_tokens)} tumour patches, {len(tumour_slides)} tumour slides", flush=True)
    self_donor = sorted(set(all_donors.patient) & set(patient_ids))
    print(f"[overlap] {len(self_donor)} normal donors are also cohort patients", flush=True)

    output = Path(args.output_dir)
    manifests: dict[str, dict] = {}
    common = dict(patient_ids=patient_ids, cancers=cancers, bags=bags, levels=levels, seed=args.seed)

    # TCGA normal slides come from the same patients as the tumour slides, so a
    # patient could otherwise be diluted with its OWN adjacent normal tissue --
    # which carries that patient's molecular state and would understate the loss.
    def exclude_self(frame: pd.DataFrame, patient: str) -> pd.DataFrame:
        return frame[frame.patient != patient] if patient in self_donor_set else frame

    self_donor_set = set(self_donor)

    if "pooled" in arms:
        manifests["pooled"] = build_arm(
            name="pooled", split=split, pool=pool, donors=donors,
            donor_lookup=lambda patient, cancer: exclude_self(donors, patient),
            output=output / "dilution_pooled.npz",
            config_extra={"donor_restriction": "none_pan_cancer_normal_pool",
                          "donor_tissue": "solid_tissue_normal",
                          "preparation_confounded": True,
                          "self_donors_excluded": self_donor}, **common)

    counts = donors.cancer.value_counts()
    test_cancers = set(cancers[split == "test"])
    matched_cancers = sorted({c for c in counts.index if counts[c] >= args.min_matched_donors} & test_cancers)
    by_cancer = {c: group for c, group in donors.groupby("cancer")}
    matched_split = np.where(np.isin(cancers, matched_cancers) & (split == "test"), "test", "train")
    if "matched" in arms:
        print(f"[matched] cancers = {matched_cancers}; test n = {(matched_split == 'test').sum()}", flush=True)
        manifests["matched"] = build_arm(
            name="matched", split=matched_split, pool=pool, donors=donors,
            donor_lookup=lambda patient, cancer: exclude_self(by_cancer.get(cancer, donors), patient),
            output=output / "dilution_matched.npz",
            config_extra={"donor_restriction": "same_cancer_type",
                          "donor_tissue": "solid_tissue_normal",
                          "preparation_confounded": True,
                          "matched_cancers": matched_cancers,
                          "min_matched_donors": args.min_matched_donors,
                          "note": "non-matched patients relabelled train so this arm keeps its own level-0 control"},
            **common)

    if "foreign_tumour" in arms:
        tumour_by_cancer = {c: group for c, group in tumour_slides.groupby("cancer")}

        def foreign_tumour_donors(patient: str, cancer: str) -> pd.DataFrame:
            group = tumour_by_cancer.get(cancer, tumour_slides)
            return group[group.patient != patient]

        manifests["foreign_tumour"] = build_arm(
            name="foreign_tumour", split=split, pool=tumour_tokens, donors=tumour_slides,
            donor_lookup=foreign_tumour_donors, output=output / "dilution_foreign_tumour.npz",
            config_extra={"donor_restriction": "same_cancer_different_patient",
                          "donor_tissue": "tumour_polygon_patches_from_the_same_store",
                          "preparation_confounded": False,
                          "purpose": "preparation-matched control: cost of uninformative patches with no domain shift"},
            **common)

    if "dx_normal" in arms and len(dx_donors):
        print(f"[dx normal] {len(dx_donors)} slides, {int(dx_donors.n.sum())} patches, "
              f"cancers {dx_donors.cancer.value_counts().to_dict()}", flush=True)
        manifests["dx_normal"] = build_arm(
            name="dx_normal", split=split, pool=all_pool, donors=dx_donors,
            donor_lookup=lambda patient, cancer: exclude_self(dx_donors, patient),
            output=output / "dilution_dx_normal.npz",
            config_extra={"donor_restriction": "none_pan_cancer_normal_pool",
                          "donor_tissue": "solid_tissue_normal_FFPE_diagnostic",
                          "preparation_confounded": False,
                          "caveat": "only 14 such slides exist in the source dataset; thin donor pool"},
            **common)

    (output / "dilution_build_manifest.json").write_text(json.dumps({
        "cohort_artifact": str(Path(args.cohort_artifact).resolve()),
        "embeddings_h5": str(Path(args.embeddings_h5).resolve()),
        "normal_staging": None if args.normal_staging is None else str(Path(args.normal_staging).resolve()),
        "normal_donor_slides": all_donors.drop(columns=["start"]).to_dict("records"),
        "arms": manifests}, indent=2, sort_keys=True), encoding="utf-8")
    print(f"[done] -> {output}", flush=True)


if __name__ == "__main__":
    main()
