"""Cohort-fit-free signature scoring, reimplemented and then *proved* against the frozen block.

``frozen_rna_targets.npz`` records its own scoring rule in ``metadata_json``:
``normalisation = {"method": "within_sample_gene_rank", "fit_population": "none"}`` and
``cohort_fit_free_target_scoring = true``.  Cohort-fit-free is the property that makes an
external cohort possible at all -- nothing about the score of an ALCHEMIST patient depends on
any TCGA patient, so the same numbers can be produced on a cohort the original code never saw.

But the *code* that produced the frozen block is not on this machine, so the rule has to be
reimplemented, and a reimplementation that is merely plausible is worthless here.  So this
module is written to be falsified: :func:`validate_against_frozen` runs the scorer on TCGA's
own raw RNA table and compares column by column with the frozen artifact.  Only signatures
that reproduce are allowed downstream, and they are dropped from *both* cohorts if they do
not, so the two sides always carry the identical target set.

Score of signature *S* in sample *j*:

    mean over g in S of  rank_j(g) / n_genes_j   -   0.5

where ``rank_j`` is the within-sample average rank (ties averaged) over every gene measured
in that sample.  Within-sample ranking is why the scorer is invariant to library size, to
the expression unit (RSEM vs TPM), and to any monotone per-sample transform -- which is what
lets TCGA's RSEM table and ALCHEMIST's STAR-Counts TPM be scored on one scale.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import rankdata

MINIMUM_REQUIRED_COVERAGE = 0.95


def read_gmt(path: str) -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                signatures[parts[0]] = [gene for gene in parts[2:] if gene]
    return signatures


def within_sample_gene_ranks(expression: np.ndarray) -> np.ndarray:
    """genes x samples -> the same shape, each column replaced by rank/n_genes in [0, 1]."""
    ranks = np.empty(expression.shape, dtype=np.float64)
    n_genes = expression.shape[0]
    for column in range(expression.shape[1]):
        ranks[:, column] = rankdata(expression[:, column], method="average") / n_genes
    return ranks


@dataclass(frozen=True)
class ScoredBlock:
    sample_ids: np.ndarray
    target_names: list[str]
    scores: np.ndarray            # (n_samples, n_targets)
    coverage: dict[str, float]
    dropped_for_coverage: list[str]


def score_signatures(expression: pd.DataFrame, signatures: dict[str, list[str]],
                     *, minimum_coverage: float = MINIMUM_REQUIRED_COVERAGE) -> ScoredBlock:
    """``expression`` is genes (index = symbol) x samples, already deduplicated by symbol."""
    ranks = within_sample_gene_ranks(expression.to_numpy(dtype=np.float64))
    position = {gene: index for index, gene in enumerate(expression.index)}
    names, columns, coverage, dropped = [], [], {}, []
    for name in sorted(signatures):
        genes = signatures[name]
        rows = [position[gene] for gene in genes if gene in position]
        fraction = len(rows) / len(genes) if genes else 0.0
        coverage[name] = fraction
        if fraction < minimum_coverage or not rows:
            dropped.append(name)
            continue
        names.append(name)
        columns.append(ranks[rows].mean(axis=0) - 0.5)
    scores = (np.vstack(columns).T if columns
              else np.zeros((expression.shape[1], 0), dtype=np.float64))
    return ScoredBlock(np.asarray(expression.columns, dtype=str), names,
                       scores.astype(np.float64), coverage, dropped)


def collapse_to_patients(sample_ids: np.ndarray, scores: np.ndarray,
                         patient_of) -> tuple[np.ndarray, np.ndarray]:
    """Average a patient's samples.  TCGA has 1,000+ patients with >1 aliquot."""
    frame = pd.DataFrame(scores)
    frame["patient"] = [patient_of(sample) for sample in sample_ids]
    grouped = frame.groupby("patient", sort=True).mean()
    return np.asarray(grouped.index, dtype=str), grouped.to_numpy(dtype=np.float64)


def validate_against_frozen(frozen_path: str, patient_ids: np.ndarray, target_names: list[str],
                            scores: np.ndarray, *, threshold: float = 0.999) -> dict:
    """Column-by-column Pearson r against the frozen TCGA block.  This is gate G1."""
    raw = np.load(frozen_path, allow_pickle=False)
    frozen_ids = raw["patient_ids"].astype(str)
    frozen_names = list(raw["target_names"].astype(str))
    frozen_scores = np.asarray(raw["scores"], dtype=np.float64)
    groups = dict(zip(frozen_names, raw["target_groups"].astype(str)))

    order = {identifier: index for index, identifier in enumerate(patient_ids)}
    rows = np.asarray([order.get(identifier, -1) for identifier in frozen_ids])
    shared = rows >= 0
    report = {"n_frozen_patients": len(frozen_ids), "n_matched_patients": int(shared.sum()),
              "threshold": threshold, "per_target": {}, "passed": [], "failed": []}
    for name in target_names:
        if name not in frozen_names:
            continue
        mine = scores[rows[shared], target_names.index(name)]
        theirs = frozen_scores[shared, frozen_names.index(name)]
        if np.std(mine) < 1e-12 or np.std(theirs) < 1e-12:
            correlation = float("nan")
        else:
            correlation = float(np.corrcoef(mine, theirs)[0, 1])
        entry = {"pearson_r": correlation,
                 "max_abs_diff": float(np.max(np.abs(mine - theirs))),
                 "group": groups.get(name, "?")}
        report["per_target"][name] = entry
        (report["passed"] if correlation >= threshold else report["failed"]).append(name)
    report["n_passed"] = len(report["passed"])
    report["n_failed"] = len(report["failed"])
    return report


def dump(report: dict) -> str:
    return json.dumps(report, indent=2, sort_keys=True)
