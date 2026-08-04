"""Build an immutable per-patient SNV target block from the MC3 public MAF.

Why this exists.  P4's stated goal is a promptable system spanning WSI, RNA,
CNV, SNV and proteomics.  Today the pipeline is WSI+RNA.  SNV was the closest
of the three missing modalities to being wired: ``calibra.known_covariate_control``
already grades mutation-derived *labels* (BRAF, KRAS, APC, PIK3CA, dMMR), but
those labels were produced by a session script that was never committed, and the
one persisted per-patient matrix on disk (``tcga_mc3_snv_binary.parquet``)
records its input sha256 and no variant filter at all -- so what a ``1`` means in
it is unknown.  An unknown filter is not a modality.  This module re-parses the
MAF with the filter written down.

What it emits is deliberately *not* a new artifact format.  It is exactly the
schema ``frozen_rna_targets.npz`` uses --- ``patient_ids``, ``cancers``,
``split``, ``target_names``, ``target_groups``, ``scores`` and ``metadata_json``
--- so ``v2.calibra.run_calibra`` and ``v2.calibra.known_covariate_control``
consume it without either of them being modified, and so the SNV block is read
through the same instrument, the same confound adjustment and the same floors as
the RNA block.  A modality that needed its own statistics would not be
comparable to the one already measured.

Three choices are load-bearing and are therefore refused rather than defaulted:

1.  **Patients with no MC3 record are excluded and named**, never imputed to
    wild-type.  Absence of sequencing is not absence of mutation, and the
    difference is invisible once it reaches a correlation.
2.  **Gene selection is fit on the development partition only.**  Which genes
    exist as targets cannot be allowed to depend on held-out patients.
3.  **Every real target gets a prevalence-matched random control** whose values
    are permuted across patients.  The marginal is preserved exactly and the
    patient pairing is destroyed, so a "channel" that the controls also clear is
    visible as such.
"""
from __future__ import annotations

import argparse
import collections
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

__all__ = [
    "PROTEIN_ALTERING",
    "ACCEPTED_FILTER_TAGS",
    "canonical_patient",
    "accepted_filter",
    "select_development_genes",
    "matched_random_controls",
    "build_snv_targets",
]

# Protein-altering coding consequences.  Silent, UTR, flank, intron, IGR and RNA
# classifications are dropped: they are real variants but they are not the thing
# "gene X is mutated" is normally taken to mean, and mixing them in would make a
# per-gene target mean two different things in two different genes.
PROTEIN_ALTERING = (
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site",
)

# MC3 ships its own FILTER tags rather than a single PASS/FAIL.  ``wga`` and
# ``native_wga_mix`` are the two the MC3 authors flag as retained-with-caveat;
# everything else (oxog, StrandBias, nonpreferredpair, broad_PoN_v2, ...) is
# dropped.  The strict-PASS-only count is recorded in the manifest so the cost of
# this choice is visible rather than argued.
ACCEPTED_FILTER_TAGS = frozenset({"PASS", "wga", "native_wga_mix"})

_MAF_COLUMNS = ("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "FILTER")


def _digest_strings(values) -> str:
    return sha256("\n".join(map(str, values)).encode("utf-8")).hexdigest()


def _digest_array(values: np.ndarray) -> str:
    return sha256(np.ascontiguousarray(values).view(np.uint8)).hexdigest()


def _file_sha256(path: str | Path, chunk_bytes: int = 4 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_patient(barcode: object) -> str:
    """Reduce a TCGA aliquot/sample barcode to its 3-field patient identifier.

    Raises on anything that is not a TCGA barcode rather than returning a
    plausible-looking string: a silently mangled identifier joins to nothing and
    shows up downstream as missing data.
    """
    text = str(barcode).strip().upper()
    pieces = text.split("-")
    if len(pieces) < 3 or pieces[0] != "TCGA" or any(not piece for piece in pieces[:3]):
        raise ValueError(f"not a TCGA barcode: {barcode!r}")
    return "-".join(pieces[:3])


def accepted_filter(value: object) -> bool:
    """True when every FILTER tag on the row is in the accepted MC3 set."""
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return False
    tags = {tag.strip() for tag in text.split(",") if tag.strip()}
    return bool(tags) and tags <= ACCEPTED_FILTER_TAGS


def select_development_genes(mutated_by_gene: dict[str, set[int]], development_rows: np.ndarray, *,
                             min_prevalence: float, min_count: int) -> list[str]:
    """Keep genes recurrent in the DEVELOPMENT partition only.

    Test rows are not consulted.  Ties are broken by gene symbol so the returned
    order is a deterministic function of the inputs and not of dict insertion.
    """
    development = set(np.flatnonzero(np.asarray(development_rows, dtype=bool)).tolist())
    if not development:
        raise ValueError("gene selection needs at least one development patient")
    kept = []
    for gene, rows in mutated_by_gene.items():
        count = len(rows & development)
        if count >= min_count and count / len(development) >= min_prevalence:
            kept.append((count, gene))
    if not kept:
        raise ValueError("no gene met the development prevalence rule; refusing to emit an empty block")
    kept.sort(key=lambda item: (-item[0], item[1]))
    return [gene for _, gene in kept]


def matched_random_controls(scores: np.ndarray, target_names, *, seed: int) -> tuple[np.ndarray, list[str]]:
    """One permuted control per real target: same marginal, no patient pairing."""
    values = np.asarray(scores, dtype=np.float32)
    rng = np.random.default_rng(seed)
    control = np.empty_like(values)
    for column in range(values.shape[1]):
        control[:, column] = values[rng.permutation(values.shape[0]), column]
    return control, [f"RANDOM_CONTROL__{name}__0" for name in target_names]


def _read_reference(reference: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Patient/cancer/split triple, taken from an existing target artifact.

    Reading the cohort off ``frozen_rna_targets.npz`` rather than re-deriving it
    guarantees the SNV block indexes the same patients, with the same cancer
    labels and the same partition assignment, as the RNA block it will be
    compared against.  Re-deriving is how two artifacts silently stop being
    comparable.
    """
    raw = np.load(Path(reference), allow_pickle=True)
    missing = {"patient_ids", "cancers", "split"} - set(raw.files)
    if missing:
        raise ValueError(f"reference artifact is missing {sorted(missing)}")
    patients = np.asarray([str(p) for p in raw["patient_ids"]])
    if len(set(patients)) != len(patients):
        raise ValueError("reference artifact has duplicate patient IDs")
    return (patients, np.asarray([str(c) for c in raw["cancers"]]),
            np.asarray([str(s) for s in raw["split"]]))


def _scan_maf(maf_path: str | Path, cohort_index: dict[str, int], *, chunk_rows: int = 500_000) -> dict:
    """Stream the MAF once, accumulating only what the targets need."""
    import pandas as pd

    protein_altering = frozenset(PROTEIN_ALTERING)
    mutated_by_gene: dict[str, set[int]] = collections.defaultdict(set)
    burden_protein = np.zeros(len(cohort_index), dtype=np.int64)
    burden_all = np.zeros(len(cohort_index), dtype=np.int64)
    aliquots_by_patient: dict[str, set[str]] = collections.defaultdict(set)
    tally = collections.Counter()

    reader = pd.read_csv(maf_path, sep="\t", usecols=list(_MAF_COLUMNS), dtype=str,
                         chunksize=chunk_rows, low_memory=False, compression="infer")
    for chunk in reader:
        tally["rows_total"] += len(chunk)
        barcodes = chunk["Tumor_Sample_Barcode"].astype(str)
        patients = barcodes.map(lambda value: "-".join(str(value).strip().upper().split("-")[:3]))
        rows = patients.map(cohort_index.get)
        in_cohort = rows.notna()
        tally["rows_in_cohort"] = tally["rows_in_cohort"] + int(in_cohort.sum())
        if not in_cohort.any():
            continue
        chunk = chunk.loc[in_cohort]
        rows = rows.loc[in_cohort].astype(np.int64).to_numpy()
        for patient, barcode in zip(patients.loc[in_cohort], barcodes.loc[in_cohort]):
            aliquots_by_patient[patient].add(str(barcode))
        keep = chunk["FILTER"].map(accepted_filter).to_numpy()
        tally["rows_filter_accepted"] = tally["rows_filter_accepted"] + int(keep.sum())
        tally["rows_filter_strict_pass"] = tally["rows_filter_strict_pass"] + int(
            (chunk["FILTER"].astype(str).str.strip() == "PASS").sum())
        if not keep.any():
            continue
        rows = rows[keep]
        classification = chunk["Variant_Classification"].astype(str).to_numpy()[keep]
        symbol = chunk["Hugo_Symbol"].astype(str).to_numpy()[keep]
        np.add.at(burden_all, rows, 1)
        coding = np.asarray([value in protein_altering for value in classification], dtype=bool)
        tally["rows_protein_altering"] = tally["rows_protein_altering"] + int(coding.sum())
        for gene, row in zip(symbol[coding], rows[coding]):
            if gene and gene != "Unknown":
                mutated_by_gene[gene].add(int(row))
        np.add.at(burden_protein, rows[coding], 1)
    return {"mutated_by_gene": dict(mutated_by_gene), "burden_protein": burden_protein,
            "burden_all": burden_all, "aliquots_by_patient": dict(aliquots_by_patient),
            "tally": dict(tally)}


def build_snv_targets(maf_path: str | Path, reference_artifact: str | Path, output_path: str | Path, *,
                      min_prevalence: float = 0.02, min_development_count: int = 100,
                      seed: int = 42, chunk_rows: int = 500_000) -> dict:
    """Parse the MC3 MAF and write ``snv_targets`` in the frozen-target schema."""
    patients, cancers, split = _read_reference(reference_artifact)
    cohort_index = {patient: row for row, patient in enumerate(patients)}
    scan = _scan_maf(maf_path, cohort_index, chunk_rows=chunk_rows)

    sequenced = np.zeros(len(patients), dtype=bool)
    for patient in scan["aliquots_by_patient"]:
        sequenced[cohort_index[patient]] = True
    if not sequenced.any():
        raise ValueError("no cohort patient appears in the MAF; check the identifier convention")
    excluded = [str(p) for p in patients[~sequenced]]

    keep_rows = np.flatnonzero(sequenced)
    remap = {int(old): new for new, old in enumerate(keep_rows)}
    mutated_by_gene = {gene: {remap[row] for row in rows if row in remap}
                       for gene, rows in scan["mutated_by_gene"].items()}
    mutated_by_gene = {gene: rows for gene, rows in mutated_by_gene.items() if rows}

    kept_patients, kept_cancers, kept_split = patients[keep_rows], cancers[keep_rows], split[keep_rows]
    development = kept_split != "test"
    genes = select_development_genes(mutated_by_gene, development,
                                     min_prevalence=min_prevalence, min_count=min_development_count)

    gene_scores = np.zeros((len(kept_patients), len(genes)), dtype=np.float32)
    for column, gene in enumerate(genes):
        gene_scores[sorted(mutated_by_gene[gene]), column] = 1.0
    burden = np.column_stack([
        np.log1p(scan["burden_protein"][keep_rows].astype(np.float64)),
        np.log1p(scan["burden_all"][keep_rows].astype(np.float64)),
    ]).astype(np.float32)

    real_names = [f"SNV_{gene}" for gene in genes] + ["SNV_BURDEN_LOG1P", "SNV_BURDEN_LOG1P_ALL"]
    real_groups = ["snv_gene"] * len(genes) + ["snv_burden"] * 2
    real_scores = np.column_stack([gene_scores, burden]).astype(np.float32)
    control_scores, control_names = matched_random_controls(real_scores, real_names, seed=seed)

    target_names = np.asarray(real_names + control_names)
    target_groups = np.asarray(real_groups + ["random_control"] * len(control_names))
    scores = np.column_stack([real_scores, control_scores]).astype(np.float32)
    if not np.isfinite(scores).all():
        raise ValueError("emitted SNV scores contain non-finite values")

    duplicates = {patient: sorted(aliquots) for patient, aliquots in scan["aliquots_by_patient"].items()
                  if len(aliquots) > 1}
    prevalence = gene_scores.mean(axis=0)
    development_prevalence = gene_scores[development].mean(axis=0)
    test_prevalence = gene_scores[~development].mean(axis=0)

    from .provenance import source_manifest

    manifest = {
        "artifact": "snv_targets",
        "schema": ["patient_ids", "cancers", "split", "target_names", "target_groups", "scores",
                   "metadata_json"],
        "source": source_manifest(configuration={
            "min_prevalence": min_prevalence, "min_development_count": min_development_count,
            "seed": seed, "protein_altering": list(PROTEIN_ALTERING),
            "accepted_filter_tags": sorted(ACCEPTED_FILTER_TAGS)}),
        "maf_path": str(Path(maf_path).resolve()), "maf_sha256": _file_sha256(maf_path),
        "reference_artifact": str(Path(reference_artifact).resolve()),
        "reference_artifact_sha256": _file_sha256(reference_artifact),
        "variant_filter": {"protein_altering": list(PROTEIN_ALTERING),
                           "accepted_filter_tags": sorted(ACCEPTED_FILTER_TAGS),
                           "maf_row_tally": scan["tally"]},
        "gene_selection": {"partition": "train+val (development only; test never consulted)",
                           "min_prevalence": float(min_prevalence),
                           "min_development_count": int(min_development_count),
                           "n_candidate_genes": len(mutated_by_gene), "n_selected_genes": len(genes),
                           "genes": genes},
        "cohort": {"reference_patient_count": int(len(patients)),
                   "emitted_patient_count": int(len(kept_patients)),
                   "excluded_no_mc3_record_count": len(excluded),
                   "excluded_no_mc3_record_patients": excluded,
                   "imputation_policy": "none; patients without an MC3 record are excluded, never set to wild-type",
                   "split_counts": {str(k): int(v) for k, v in
                                    zip(*np.unique(kept_split, return_counts=True))},
                   "cancer_count": int(len(set(map(str, kept_cancers)))),
                   "patients_with_multiple_aliquots": len(duplicates),
                   "multiple_aliquot_examples": dict(sorted(duplicates.items())[:20])},
        "targets": {"n_real": len(real_names), "n_random_control": len(control_names),
                    "groups": {str(k): int(v) for k, v in
                               zip(*np.unique(target_groups, return_counts=True))},
                    "gene_prevalence_overall": {gene: float(v) for gene, v in zip(genes, prevalence)},
                    "gene_prevalence_development": {gene: float(v) for gene, v in
                                                    zip(genes, development_prevalence)},
                    "gene_prevalence_test": {gene: float(v) for gene, v in zip(genes, test_prevalence)},
                    "burden_protein_altering_mean": float(scan["burden_protein"][keep_rows].mean()),
                    "burden_all_mean": float(scan["burden_all"][keep_rows].mean())},
        "digests": {"patient_ids": _digest_strings(kept_patients), "cancers": _digest_strings(kept_cancers),
                    "split": _digest_strings(kept_split), "target_names": _digest_strings(target_names),
                    "target_groups": _digest_strings(target_groups), "scores": _digest_array(scores)},
    }

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, patient_ids=kept_patients.astype(str), cancers=kept_cancers.astype(str),
                        split=kept_split.astype(str), target_names=target_names.astype(str),
                        target_groups=target_groups.astype(str), scores=scores,
                        metadata_json=np.asarray(json.dumps(manifest, sort_keys=True)))
    manifest["output_path"] = str(output.resolve())
    manifest["output_sha256"] = _file_sha256(output)
    output.with_suffix(output.suffix + ".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maf", required=True, help="MC3 public MAF (.maf.gz)")
    parser.add_argument("--reference-artifact", required=True,
                        help="existing target artifact supplying patient_ids/cancers/split")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-prevalence", type=float, default=0.02)
    parser.add_argument("--min-development-count", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chunk-rows", type=int, default=500_000)
    args = parser.parse_args()
    manifest = build_snv_targets(args.maf, args.reference_artifact, args.output,
                                 min_prevalence=args.min_prevalence,
                                 min_development_count=args.min_development_count,
                                 seed=args.seed, chunk_rows=args.chunk_rows)
    print(json.dumps({k: manifest[k] for k in ("output_path", "output_sha256")}, indent=2), flush=True)
    print(f"[snv] patients={manifest['cohort']['emitted_patient_count']} "
          f"excluded={manifest['cohort']['excluded_no_mc3_record_count']} "
          f"genes={manifest['gene_selection']['n_selected_genes']} "
          f"targets={manifest['targets']['n_real']}+{manifest['targets']['n_random_control']} controls",
          flush=True)


if __name__ == "__main__":
    main()
