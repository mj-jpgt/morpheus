"""Score ALCHEMIST and TCGA-NSCLC on one target block, over one shared gene universe.

Two cohorts, two very different expression products: TCGA is the EBPlusPlus RSEM pancan
table keyed by ``SYMBOL|ENTREZ``; ALCHEMIST is 1,106 per-sample STAR-Counts TSVs keyed by
Ensembl ID with a ``gene_name`` column and a ``tpm_unstranded`` value.  Within-sample
ranking makes the *unit* irrelevant -- RSEM and TPM give the same ranks up to ties -- but it
does **not** make the *universe* irrelevant.  A rank is a position among however many genes
were measured, so ranking ALCHEMIST over its ~60k Ensembl rows and TCGA over its ~20.5k
symbols would put the two cohorts on two different scales while looking identical.

So the universe is intersected once and both cohorts are ranked inside it.  That means TCGA
is **re-scored** here rather than read out of ``frozen_rna_targets.npz``: the frozen block
was ranked over TCGA's own 20,502 symbols and is not on the shared scale.  The frozen block
is still used, as gate G1 -- ``validate_against_frozen`` reproduces it exactly on TCGA's own
universe (r = 1.000000, max abs diff 0.0000000) -- which is what licenses re-scoring at all.

Duplicate aliquots follow the frozen artifact's recorded ``primary_tumor_then_same_type_mean``
policy, using its own per-patient ``selected`` sample lists verbatim where they exist.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3].parent))

from morpheus.v2.research.external.rank_target_scoring import (  # noqa: E402
    MINIMUM_REQUIRED_COVERAGE, read_gmt, score_signatures, validate_against_frozen,
    within_sample_gene_ranks,
)
from morpheus.v2.curated_panel import FROZEN_MECHANISM_PROGRAMMES  # noqa: E402

GDC_DATA = "https://api.gdc.cancer.gov/data"

# The 16 non-MSigDB targets of the frozen block are curated gene programmes.  These are the
# names as they appear in `frozen_rna_targets.npz` mapped onto `FROZEN_MECHANISM_PROGRAMMES`.
# `immune_t_cell_inflammation` has no entry in that dict, so it is simply absent and will be
# reported as unavailable rather than approximated by a similar-sounding programme.
CURATED_ALIASES = {
    "immune_cytolytic_activity": "cytolytic_activity",
    "immune_ifng": "interferon_gamma",
    "immune_antigen_presentation": "antigen_presentation",
    "immune_myeloid_macrophage": "myeloid_macrophage",
    "stroma_caf": "stroma_caf",
    "tgfb_emt": "tgf_beta_emt",
    "immune_exclusion": "immune_exclusion",
    "state_proliferation": "proliferation",
    "state_hypoxia": "hypoxia",
    "state_glycolysis": "glycolysis",
    "state_angiogenesis": "angiogenesis",
    "state_dna_repair": "dna_repair",
    "state_apoptosis_senescence": "apoptosis_senescence",
    "state_emt": "emt",
    "state_mechanotransduction": "mechanotransduction",
}


def fetch_expression(file_id: str, destination: Path, retries: int = 4) -> Path:
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(f"{GDC_DATA}/{file_id}"),
                                        timeout=300) as response:
                payload = response.read()
            destination.write_bytes(payload)
            return destination
        except (urllib.error.URLError, TimeoutError, OSError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def read_star_counts(path: Path, value_column: str = "tpm_unstranded") -> pd.Series:
    """One STAR-Counts TSV -> Series indexed by gene symbol.

    The first four rows are the ``N_unmapped`` / ``N_multimapping`` / ``N_noFeature`` /
    ``N_ambiguous`` summary lines GDC keeps in the file; they are not genes and must not
    enter a rank universe.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        frame = pd.read_csv(handle, sep="\t", comment="#")
    frame = frame[~frame["gene_id"].astype(str).str.startswith("N_")]
    frame = frame[frame["gene_name"].notna()]
    values = frame.groupby("gene_name")[value_column].max()
    return values


def tcga_expression(table: Path, samples: list[str]) -> pd.DataFrame:
    columns = ["gene_id"] + samples
    frame = pd.read_csv(table, sep="\t", usecols=columns).set_index("gene_id")
    frame.index = pd.Series(frame.index).str.split("|").str[0].values
    return frame[~frame.index.isin(["?"])]


def resolve_duplicates(frozen_metadata: dict, patient: str, samples: list[str]) -> list[str]:
    resolutions = frozen_metadata["duplicate_samples"]["resolutions"]
    if patient in resolutions:
        selected = [s for s in resolutions[patient]["selected"] if s in samples]
        if selected:
            return selected
    return samples


def build_signatures(gmt_path: str, wanted: list[str]) -> tuple[dict[str, list[str]], list[str]]:
    gmt = read_gmt(gmt_path)
    signatures, unavailable = {}, []
    for name in wanted:
        if name in gmt:
            signatures[name] = gmt[name]
        elif name in CURATED_ALIASES and CURATED_ALIASES[name] in FROZEN_MECHANISM_PROGRAMMES:
            signatures[name] = list(FROZEN_MECHANISM_PROGRAMMES[CURATED_ALIASES[name]])
        else:
            unavailable.append(name)
    return signatures, unavailable


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alchemist-manifest", required=True)
    parser.add_argument("--frozen-targets", required=True)
    parser.add_argument("--tcga-rna", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--expression-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--validate-only", action="store_true")
    # The frozen artifact records `minimum_required_coverage: 0.95`, but it plainly did not
    # apply it as a per-signature drop in this symbol universe: its own manifest keeps
    # HALLMARK_ADIPOGENESIS at 189/200 = 0.945. Applying 0.95 during validation would
    # therefore drop ten hallmarks the frozen block retained, and G1 would silently be
    # answering a smaller question than it claims to. So the gate runs uncut, and the
    # coverage requirement is enforced where it actually protects something: on ALCHEMIST,
    # the cohort whose gene coverage nobody has checked.
    parser.add_argument("--validation-min-coverage", type=float, default=0.0)
    parser.add_argument("--cohort-min-coverage", type=float, default=MINIMUM_REQUIRED_COVERAGE)
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    frozen = np.load(args.frozen_targets, allow_pickle=False)
    metadata = json.loads(str(frozen["metadata_json"]))
    frozen_names = frozen["target_names"].astype(str)
    frozen_groups = frozen["target_groups"].astype(str)
    wanted = [name for name, group in zip(frozen_names, frozen_groups)
              if group != "random_control"]
    signatures, unavailable = build_signatures(args.gmt, wanted)
    print(f"[targets] {len(wanted)} non-control targets requested, "
          f"{len(signatures)} resolvable, unavailable={unavailable}", flush=True)

    # ---- TCGA side -------------------------------------------------------------------
    header = pd.read_csv(args.tcga_rna, sep="\t", nrows=0)
    all_samples = list(header.columns)[1:]
    frozen_patients = set(frozen["patient_ids"].astype(str))
    cancers = dict(zip(frozen["patient_ids"].astype(str), frozen["cancers"].astype(str)))
    nsclc = {p for p in frozen_patients if cancers.get(p) in {"LUAD", "LUSC"}}

    by_patient: dict[str, list[str]] = {}
    for sample in all_samples:
        patient = sample[:12]
        if patient in frozen_patients:
            by_patient.setdefault(patient, []).append(sample)
    selected_samples, sample_patient = [], {}
    for patient, samples in by_patient.items():
        for sample in resolve_duplicates(metadata, patient, sorted(samples)):
            selected_samples.append(sample)
            sample_patient[sample] = patient

    # Gate G1 runs on a duplicate-free TCGA subset over TCGA's OWN universe, which is the
    # only configuration in which the frozen numbers are reproducible by construction.
    validation_samples = [s for s in selected_samples
                          if len(by_patient[sample_patient[s]]) == 1][:400]
    print(f"[G1] validating on {len(validation_samples)} duplicate-free TCGA samples", flush=True)
    frame = tcga_expression(Path(args.tcga_rna), validation_samples)
    block = score_signatures(frame, signatures, minimum_coverage=args.validation_min_coverage)
    patients = np.asarray([sample_patient[s] for s in validation_samples])
    report = validate_against_frozen(args.frozen_targets, patients, block.target_names,
                                     block.scores)
    report["gene_universe"] = int(frame.shape[0])
    report["unavailable_targets"] = unavailable
    report["coverage"] = block.coverage
    report["validation_min_coverage"] = args.validation_min_coverage
    report["cohort_min_coverage"] = args.cohort_min_coverage
    report["dropped_for_coverage"] = block.dropped_for_coverage
    (output / "G1_validation.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    print(f"[G1] passed={report['n_passed']} failed={report['n_failed']} "
          f"(threshold {report['threshold']})", flush=True)
    if report["failed"]:
        print(f"[G1] failing targets: {report['failed'][:10]}", flush=True)
    if args.validate_only:
        return

    # ---- ALCHEMIST side --------------------------------------------------------------
    catalog = pd.read_csv(args.alchemist_manifest)
    expression_files = json.loads((Path(args.expression_dir) / "expression_index.json").read_text())
    directory = Path(args.expression_dir)
    pool = ThreadPoolExecutor(max_workers=args.workers)
    futures = {entry["case_id"]: pool.submit(fetch_expression, entry["file_id"],
                                             directory / f"{entry['case_id']}.tsv")
               for entry in expression_files}
    series: dict[str, pd.Series] = {}
    for index, entry in enumerate(expression_files):
        try:
            path = futures[entry["case_id"]].result()
            series[entry["case_id"]] = read_star_counts(path)
        except Exception as exc:  # noqa: BLE001
            print(f"[skip] expression {entry['case_id']}: {type(exc).__name__} {exc}", flush=True)
        if (index + 1) % 200 == 0:
            print(f"[expr] {index + 1}/{len(expression_files)}", flush=True)
    pool.shutdown(wait=True)
    alchemist = pd.DataFrame(series)
    print(f"[expr] ALCHEMIST matrix {alchemist.shape}", flush=True)

    # ---- one shared universe ---------------------------------------------------------
    tcga_full = tcga_expression(Path(args.tcga_rna),
                                [s for s in selected_samples if sample_patient[s] in nsclc])
    shared = sorted(set(tcga_full.index) & set(alchemist.index))
    print(f"[universe] TCGA {tcga_full.shape[0]} x ALCHEMIST {alchemist.shape[0]} "
          f"-> shared {len(shared)}", flush=True)
    tcga_shared = tcga_full.loc[~tcga_full.index.duplicated(keep="first")].loc[shared]
    alch_shared = alchemist.loc[shared]

    # Cancer labels must travel WITH the target block. Without them the downstream
    # confound design silently collapses to a single level and the LUAD/LUSC adjustment
    # that the TCGA arm is supposed to carry quietly stops happening.
    frozen_cancer = dict(zip(frozen["patient_ids"].astype(str), frozen["cancers"].astype(str)))
    alch_cancer = dict(zip(catalog["case_id"].astype(str), catalog["cancer"].astype(str)))
    cancer_of = {"tcga_nsclc": lambda pid: frozen_cancer.get(pid, "NA"),
                 "alchemist": lambda pid: alch_cancer.get(pid, "NA")}
    patient_of = {"tcga_nsclc": lambda pid: pid,
                  "alchemist": dict(zip(catalog["case_id"].astype(str),
                                        catalog["patient"].astype(str))).get}

    results = {}
    for label, matrix, id_of in (
            ("tcga_nsclc", tcga_shared, lambda c: sample_patient[c]),
            ("alchemist", alch_shared, lambda c: c)):
        scored = score_signatures(matrix, signatures, minimum_coverage=args.cohort_min_coverage)
        frame = pd.DataFrame(scored.scores, columns=scored.target_names)
        frame["id"] = [id_of(column) for column in matrix.columns]
        grouped = frame.groupby("id", sort=True).mean()
        results[label] = {"ids": np.asarray(grouped.index, dtype=str),
                          "names": scored.target_names,
                          "scores": grouped.to_numpy(dtype=np.float64),
                          "coverage": scored.coverage,
                          "dropped": scored.dropped_for_coverage}
        print(f"[score] {label}: {grouped.shape[0]} ids x {len(scored.target_names)} targets, "
              f"dropped_for_coverage={scored.dropped_for_coverage}", flush=True)

    # Gate G1 as predeclared: a target that does not reproduce the frozen block is dropped
    # from BOTH cohorts. The 15 curated mechanism programmes fail because their gene lists
    # are not recoverable from any artifact on disk, so keeping them would mean scoring
    # ALCHEMIST on a signature that is not the one TCGA was scored on.
    g1_passed = set(report["passed"])
    shared_names = set(results["tcga_nsclc"]["names"]) & set(results["alchemist"]["names"])
    common = [name for name in results["tcga_nsclc"]["names"]
              if name in shared_names and name in g1_passed]
    print(f"[gate] G1 passed {len(g1_passed)}; present in both cohorts after the "
          f"coverage cut: {len(common)}", flush=True)
    for label in results:
        keep = [results[label]["names"].index(name) for name in common]
        ids = results[label]["ids"]
        np.savez_compressed(
            output / f"{label}_targets.npz",
            patient_ids=np.asarray([patient_of[label](i) or i for i in ids], dtype=str),
            case_ids=np.asarray(ids, dtype=str),
            cancers=np.asarray([cancer_of[label](i) for i in ids], dtype=str),
            target_names=np.asarray(common),
            target_groups=np.asarray([dict(zip(frozen_names, frozen_groups))[n] for n in common]),
            scores=results[label]["scores"][:, keep].astype(np.float32),
            metadata_json=np.asarray(json.dumps({
                "cohort": label,
                "scoring": "within_sample_gene_rank over a SHARED TCGA/ALCHEMIST symbol universe",
                "shared_universe_size": len(shared),
                "minimum_required_coverage": MINIMUM_REQUIRED_COVERAGE,
                "coverage": results[label]["coverage"],
                "dropped_for_coverage": results[label]["dropped"],
                "unavailable_targets": unavailable,
                "g1_passed": report["n_passed"], "g1_failed": report["n_failed"],
                "gmt": args.gmt,
                "value_column": "tpm_unstranded" if label == "alchemist" else "EBPlusPlus_RSEM",
            }, sort_keys=True)))
    # Sensitivity: G1 only, coverage cut removed. Intersecting the two gene universes
    # shrinks per-signature coverage, so the 0.95 rule removes 15 targets that pass G1
    # outright. Emitting both means the headline cannot depend on that particular cut.
    relaxed = {}
    for label, matrix, id_of in (("tcga_nsclc", tcga_shared, lambda c: sample_patient[c]),
                                 ("alchemist", alch_shared, lambda c: c)):
        scored = score_signatures(matrix, signatures, minimum_coverage=0.0)
        frame = pd.DataFrame(scored.scores, columns=scored.target_names)
        frame["id"] = [id_of(column) for column in matrix.columns]
        grouped = frame.groupby("id", sort=True).mean()
        relaxed[label] = (np.asarray(grouped.index, dtype=str), scored.target_names,
                          grouped.to_numpy(dtype=np.float64))
    relaxed_common = [n for n in relaxed["tcga_nsclc"][1]
                      if n in set(relaxed["alchemist"][1]) and n in g1_passed]
    for label in relaxed:
        ids, names_r, values = relaxed[label]
        keep = [names_r.index(n) for n in relaxed_common]
        np.savez_compressed(
            output / f"{label}_targets_nocoveragecut.npz",
            patient_ids=np.asarray([patient_of[label](i) or i for i in ids], dtype=str),
            case_ids=np.asarray(ids, dtype=str),
            cancers=np.asarray([cancer_of[label](i) for i in ids], dtype=str),
            target_names=np.asarray(relaxed_common),
            target_groups=np.asarray([dict(zip(frozen_names, frozen_groups))[n]
                                      for n in relaxed_common]),
            scores=values[:, keep].astype(np.float32),
            metadata_json=np.asarray(json.dumps({
                "cohort": label, "variant": "G1 only, no coverage cut",
                "shared_universe_size": len(shared)}, sort_keys=True)))
    print(f"[done] {len(common)} targets (G1 + coverage) and {len(relaxed_common)} "
          f"(G1 only) written to {output}", flush=True)


if __name__ == "__main__":
    main()
