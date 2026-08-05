"""TEST 2 -- do P3's 29 certified causal NAMES predict which untrained targets an axis reads?

Predeclared in
``NOTEBOOK_ENTRIES/PREDECLARED_p4_composed_readout_and_causal_name_bridge_20260805T0650Z.md``,
committed (``b83ac7c``) before this file was written. The process-family map of §2.3 is
transcribed here from that file verbatim; it was fixed before any ``r(j, t)`` existed.

Hypothesis. An axis named for a biological process should read queries about that process and
about processes adjacent to it, including processes never trained on. If names carry no
generalisation, name-adjacency will not predict readout strength once the axis and target main
effects are removed.

**No statistic is defined here.** Every number comes from ``v2/calibra`` and ``v2``:

* ``spectral.heldout_cca_projection``   -- the image direction FITTED for PCA axis j, on train rows
* ``spectral.paired_absolute_correlation`` -- that projection read against target t, on test rows
* ``residualise.confound_design`` / ``cross_fitted_residuals`` / ``pooled_tissue_source_site``
* ``perturbation_basis_common.load_aligned_response`` -- the K562 atom responses, ``tcga_sd``
  scaled exactly as ``causal_attribution.attribution_report`` scales them
* ``causal_attribution.atom_cosines``   -- the name-vs-gene-set adjacency
* ``baseline_target_common.load_reference_targets`` / ``development_expression_moments``

The two nuisances, and how they are removed
-------------------------------------------
Legible axes read everything and some targets are readable by everything. The 29 x 40 matrix of
``r(j, t)`` is **double-centred** (row means, then column means, then grand mean added back)
before any association is taken, and the permutation null shuffles adjacency **across the 40
targets within each axis**, which preserves both main effects exactly and destroys only the
name-target correspondence.

The matched-legibility control
------------------------------
Certified axes are compared with uncertified axes of similar legibility -- 1:1 nearest neighbour
on ``legibility__d2_h_seed42__wsi_biology``, without replacement, caliper 0.02, greedy in
descending certified legibility -- or the test would merely rediscover that legible axes read
things well.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.baseline_target_common import (development_expression_moments,
                                                load_reference_targets)
from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import heldout_cca_projection, paired_absolute_correlation
from morpheus.v2.causal_attribution import atom_cosines
from morpheus.v2.curated_panel import FROZEN_MECHANISM_PROGRAMMES
from morpheus.v2.perturbation_basis_common import load_aligned_response

LEGIBILITY = "legibility__d2_h_seed42__wsi_biology"
UNTRAINED_GROUPS = ("heldout_pathway", "immune_tme", "tumour_state")
CALIPER = 0.02
N_COMPONENTS = 32

#: §2.3 of the predeclaration, transcribed verbatim. Assigned from ``top_perturbed_genes``
#: alone, before any readout statistic existed.
AXIS_FAMILY = {
    "PCA_004": "RIBO", "PCA_024": "RIBO", "PCA_075": "RIBO", "PCA_078": "RIBO",
    "PCA_013": "MITO", "PCA_094": "MITO", "PCA_097": "MITO", "PCA_098": "MITO",
    "PCA_102": "MITO",
    "PCA_031": "REPL", "PCA_032": "REPL", "PCA_085": "REPL", "PCA_108": "REPL",
    "PCA_111": "REPL",
    "PCA_038": "TXN", "PCA_068": "TXN", "PCA_104": "TXN", "PCA_109": "TXN", "PCA_112": "TXN",
    "PCA_122": "TXN", "PCA_126": "TXN",
    "PCA_047": "SECR", "PCA_051": "SECR", "PCA_063": "SECR", "PCA_072": "SECR",
    "PCA_080": "SPLC", "PCA_114": "SPLC",
    "PCA_007": "OTHER", "PCA_027": "OTHER",
}

#: §2.3 of the predeclaration, transcribed verbatim. Assigned from target names alone.
TARGET_FAMILY = {
    "KEGG_MEDICUS_PATHOGEN_ARSENIC_TO_ELECTRON_TRANSFER_IN_COMPLEX_II": "MITO",
    "KEGG_MEDICUS_ENV_FACTOR_ARSENIC_TO_ELECTRON_TRANSFER_IN_COMPLEX_IV": "MITO",
    "state_glycolysis": "MITO", "state_hypoxia": "MITO",
    "KEGG_MEDICUS_ENV_FACTOR_NNK_NNN_TO_CHRNA7_E2F_SIGNALING_PATHWAY": "REPL",
    "KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS": "REPL",
    "KEGG_MEDICUS_ENV_FACTOR_DCE_TO_DNA_ADDUCTS": "REPL",
    "state_proliferation": "REPL", "state_dna_repair": "REPL",
    "KEGG_MEDICUS_PATHOGEN_EBV_EBNA1_TO_P53_MEDIATED_TRANSCRIPTION": "TXN",
    "KEGG_MEDICUS_PATHOGEN_EBV_EBNA2_TO_RBP_JK_MEDIATED_TRANSCRIPTION": "TXN",
    "KEGG_MEDICUS_ENV_FACTOR_TCDD_TO_AHR_SIGNALING_PATHWAY": "TXN",
    "KEGG_MEDICUS_ENV_FACTOR_E2_TO_NUCLEAR_INITIATED_ESTROGEN_SIGNALING_PATHWAY": "TXN",
    "KEGG_MEDICUS_PATHOGEN_EBV_BARF1_TO_INTRINSIC_APOPTOTIC_PATHWAY": "APOP",
    "KEGG_MEDICUS_ENV_FACTOR_PARAQUAT_TO_FAS_JNK_SIGNALING_PATHWAY": "APOP",
    "state_apoptosis_senescence": "APOP",
    "KEGG_MEDICUS_ENV_FACTOR_BENZO_A_PYRENRE_TO_CYP_MEDIATED_METABOLISM": "XENO",
    "KEGG_MEDICUS_ENV_FACTOR_METALS_TO_KEAP1_NRF2_SIGNALIG_PATHWAY": "XENO",
    "KEGG_MEDICUS_ENV_FACTOR_E2_TO_RAS_ERK_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_METALS_TO_RAS_ERK_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_NNK_NNN_TO_RAS_ERK_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_METALS_TO_JNK_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_METALS_TO_NFKB_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_NICOTINE_NNK_TO_PI3K_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_NNK_NNN_TO_PI3K_SIGNALING_PATHWAY_N01339": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_NNK_NNN_TO_PI3K_SIGNALING_PATHWAY_N01350": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_NICOTINE_TO_JAK_STAT_SIGNALING_PATHWAY": "SIGNAL",
    "KEGG_MEDICUS_ENV_FACTOR_IRON_TO_ANTEROGRADE_AXONAL_TRANSPORT": "TRANSPORT",
    "KEGG_MEDICUS_ENV_FACTOR_ZN_TO_ANTEROGRADE_AXONAL_TRANSPORT": "TRANSPORT",
    "state_angiogenesis": "OTHER_STATE", "state_emt": "OTHER_STATE",
    "state_mechanotransduction": "OTHER_STATE",
}
#: Every ``immune_tme`` target maps to IMMUNE; filled in at load time from the group column.

#: §2.3, verbatim: the three explicitly adjacent cross-family pairs. Everything else is
#: non-adjacent, including every SIGNAL / IMMUNE / TRANSPORT / XENO / OTHER_STATE target,
#: which have no adjacent certified axis family at all.
ADJACENT_PAIRS = {frozenset({"RIBO", "REPL"}), frozenset({"SECR", "MITO"}),
                  frozenset({"REPL", "APOP"})}

#: The 16 non-MSigDB targets are curated programmes, not GMT sets. Alias map copied from
#: ``research/external/build_alchemist_targets.CURATED_ALIASES``; ``immune_t_cell_inflammation``
#: has no curated programme and is reported as having no gene set rather than approximated.
CURATED_ALIASES = {
    "immune_cytolytic_activity": "cytolytic_activity", "immune_ifng": "interferon_gamma",
    "immune_antigen_presentation": "antigen_presentation",
    "immune_myeloid_macrophage": "myeloid_macrophage", "stroma_caf": "stroma_caf",
    "tgfb_emt": "tgf_beta_emt", "immune_exclusion": "immune_exclusion",
    "state_proliferation": "proliferation", "state_hypoxia": "hypoxia",
    "state_glycolysis": "glycolysis", "state_angiogenesis": "angiogenesis",
    "state_dna_repair": "dna_repair", "state_apoptosis_senescence": "apoptosis_senescence",
    "state_emt": "emt", "state_mechanotransduction": "mechanotransduction",
}


def read_gmt(path: str) -> dict:
    """One reader, the same three-field GMT convention as ``rank_target_scoring.read_gmt``."""
    out = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                out[parts[0]] = [gene for gene in parts[2:] if gene]
    return out


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return float(spearmanr(a[ok], b[ok]).statistic)


def _double_centre(matrix: np.ndarray) -> np.ndarray:
    m = np.asarray(matrix, dtype=np.float64)
    return m - m.mean(axis=1, keepdims=True) - m.mean(axis=0, keepdims=True) + m.mean()


def _association(readout: np.ndarray, adjacency: np.ndarray) -> float:
    return _spearman(_double_centre(readout).ravel(), _double_centre(adjacency).ravel())


def _permutation_null(readout: np.ndarray, adjacency: np.ndarray, *, n_draws: int,
                      seed: int) -> dict:
    """Shuffle adjacency ACROSS TARGETS WITHIN EACH AXIS. Both main effects survive exactly."""
    rng = np.random.default_rng(seed)
    observed = _association(readout, adjacency)
    null = np.empty(n_draws)
    for d in range(n_draws):
        permuted = np.stack([row[rng.permutation(adjacency.shape[1])] for row in adjacency])
        null[d] = _association(readout, permuted)
    exceed = int(np.sum(null >= observed))
    return {"observed": float(observed), "null_median": float(np.median(null)),
            "null_p95": float(np.percentile(null, 95)),
            "permutation_p": float((exceed + 1) / (n_draws + 1)), "n_draws": int(n_draws)}


def _match_on_legibility(table: pd.DataFrame) -> dict:
    """1:1 nearest neighbour on legibility, without replacement, caliper 0.02, greedy."""
    certified = table[table["causal_name_certified"] == True]  # noqa: E712
    pool = table[table["causal_name_certified"] != True].copy()  # noqa: E712
    order = certified.sort_values(LEGIBILITY, ascending=False)
    available = list(pool.index)
    pairs, dropped = [], []
    for index, row in order.iterrows():
        if not available:
            dropped.append(str(row["axis"]))
            continue
        distances = np.abs(pool.loc[available, LEGIBILITY].to_numpy(float) - float(row[LEGIBILITY]))
        best = int(np.argmin(distances))
        if distances[best] > CALIPER:
            dropped.append(str(row["axis"]))
            continue
        partner = available.pop(best)
        pairs.append({"certified_axis": str(row["axis"]),
                      "certified_legibility": float(row[LEGIBILITY]),
                      "control_axis": str(pool.loc[partner, "axis"]),
                      "control_legibility": float(pool.loc[partner, LEGIBILITY]),
                      "delta": float(distances[best])})
    deltas = np.asarray([p["delta"] for p in pairs], dtype=float)
    return {"pairs": pairs, "dropped_certified_axes": dropped,
            "n_pairs": len(pairs), "n_dropped": len(dropped),
            "mean_abs_delta": float(deltas.mean()) if len(deltas) else float("nan"),
            "max_abs_delta": float(deltas.max()) if len(deltas) else float("nan"),
            "balance_ok": bool(len(deltas) and deltas.mean() <= 0.005 and deltas.max() <= CALIPER)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True)
    parser.add_argument("--state", default="wsi_biology")
    parser.add_argument("--pca-targets", required=True)
    parser.add_argument("--rna-targets", required=True)
    parser.add_argument("--attribution", required=True)
    parser.add_argument("--pbs-targets", required=True)
    parser.add_argument("--rna-table", required=True)
    parser.add_argument("--perturbation", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--partition", default="test")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--split-seed", type=int, default=42)
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # ---- the image state, on ``_axis_legibility``'s own footing -------------------------
    raw = np.load(args.artifact, allow_pickle=True)
    patient_ids = np.asarray([str(p) for p in raw["patient_ids"]])
    cancers = np.asarray([str(c) for c in raw["cancers"]])
    split = np.asarray([str(s) for s in raw["split"]])
    mask = np.ones(len(split), bool) if args.partition == "all" else (split == args.partition)

    pca = np.load(args.pca_targets, allow_pickle=True)
    pca_ids = {str(p): i for i, p in enumerate(np.asarray(pca["patient_ids"]))}
    rna = np.load(args.rna_targets, allow_pickle=True)
    rna_ids = {str(p): i for i, p in enumerate(np.asarray(rna["patient_ids"]))}
    mask &= np.asarray([p in pca_ids and p in rna_ids for p in patient_ids])

    rows = np.flatnonzero(mask)
    tss, _ = pooled_tissue_source_site(patient_ids[rows])
    design = confound_design(pd.DataFrame({"cancer": cancers[rows], "tss": tss}), ["cancer", "tss"])
    features = cross_fitted_residuals(np.asarray(raw[args.state], dtype=np.float64)[rows],
                                      design, seed=args.seed)
    pca_scores = cross_fitted_residuals(
        np.asarray(pca["scores"], dtype=np.float64)[[pca_ids[p] for p in patient_ids[rows]]],
        design, seed=args.seed)

    target_names = np.asarray([str(n) for n in rna["target_names"]])
    target_groups = np.asarray([str(g) for g in rna["target_groups"]])
    keep = np.flatnonzero(np.isin(target_groups, list(UNTRAINED_GROUPS)))
    target_names, target_groups = target_names[keep], target_groups[keep]
    target_scores = cross_fitted_residuals(
        np.asarray(rna["scores"], dtype=np.float64)[
            np.ix_([rna_ids[p] for p in patient_ids[rows]], keep)], design, seed=args.seed)
    print(f"[state] n={len(rows)} axes={pca_scores.shape[1]} targets={len(target_names)}",
          flush=True)

    # ---- r(j, t): axis j's FITTED image direction, read against target t on held-out rows
    order = np.random.default_rng(args.split_seed).permutation(len(rows))
    cut = len(rows) // 2
    train, test = order[:cut], order[cut:]
    n_axes = pca_scores.shape[1]
    readout = np.full((n_axes, len(target_names)), np.nan)
    self_readout = np.full(n_axes, np.nan)
    for j in range(n_axes):
        px, py = heldout_cca_projection(features, pca_scores[:, j][:, None], train, test,
                                        n_components=N_COMPONENTS)
        if px.size == 0:
            continue
        self_readout[j] = paired_absolute_correlation(px, py)
        for t in range(len(target_names)):
            readout[j, t] = paired_absolute_correlation(px, target_scores[test, t])
    print("[readout] done", flush=True)

    # ---- adjacency A1: does perturbing the axis's named genes move the target's genes? ---
    reference = load_reference_targets(args.pbs_targets)
    _expression, _transform, _mean, scale = development_expression_moments(args.rna_table,
                                                                          reference)
    genes = [str(g) for g in reference["genes"]]
    primary = load_aligned_response(args.perturbation, genes, scale, scaling="tcga_sd")
    aligned_genes = np.asarray(genes)[primary["gene_index"]]
    gene_position = {g: i for i, g in enumerate(aligned_genes)}

    gmt = read_gmt(args.gmt)
    gene_sets, missing_sets = {}, []
    for name in target_names:
        if name in gmt:
            members = gmt[name]
        elif name in CURATED_ALIASES:
            members = list(FROZEN_MECHANISM_PROGRAMMES[CURATED_ALIASES[name]])
        else:
            missing_sets.append(name)
            continue
        hit = [gene_position[g] for g in members if g in gene_position]
        if hit:
            gene_sets[name] = np.asarray(hit, dtype=np.int64)
        else:
            missing_sets.append(name)

    directions = np.zeros((len(aligned_genes), len(target_names)))
    for t, name in enumerate(target_names):
        if name in gene_sets:
            indicator = np.zeros(len(aligned_genes))
            indicator[gene_sets[name]] = 1.0
            directions[:, t] = indicator - indicator.mean()
    cosine = np.abs(atom_cosines(primary["response"], directions))     # [atom, target]
    atom_row = {str(a): i for i, a in enumerate(primary["atom_ids"])}

    table = pd.read_csv(args.attribution)
    adjacency_a1 = np.full((n_axes, len(target_names)), np.nan)
    adjacency_a1_lit = np.full((n_axes, len(target_names)), np.nan)
    universe = len(aligned_genes)
    from scipy.stats import hypergeom
    for _, record in table.iterrows():
        j = int(record["axis_index"])
        atoms = [atom_row[a] for a in str(record["top_atoms"]).split(";") if a in atom_row]
        named = [g for g in str(record["top_perturbed_genes"]).split(";") if g in gene_position]
        for t, name in enumerate(target_names):
            if name not in gene_sets:
                continue
            if atoms:
                adjacency_a1[j, t] = float(np.mean(cosine[atoms, t]))
            members = set(gene_sets[name].tolist())
            overlap = sum(1 for g in named if gene_position[g] in members)
            adjacency_a1_lit[j, t] = float(hypergeom.sf(overlap - 1, universe, len(members),
                                                        max(len(named), 1)))

    # ---- A2: the predeclared process-family map -----------------------------------------
    family_of_target = {}
    for name, group in zip(target_names, target_groups):
        family_of_target[name] = ("IMMUNE" if group == "immune_tme"
                                  else TARGET_FAMILY.get(name, "UNMAPPED"))
    unmapped = sorted({n for n, f in family_of_target.items() if f == "UNMAPPED"})
    adjacency_a2 = np.zeros((n_axes, len(target_names)))
    axis_name = {int(r["axis_index"]): str(r["axis"]) for _, r in table.iterrows()}
    for j in range(n_axes):
        axis_family = AXIS_FAMILY.get(axis_name.get(j, ""), None)
        for t, name in enumerate(target_names):
            target_family = family_of_target[name]
            if axis_family is None:
                adjacency_a2[j, t] = np.nan
            elif axis_family == target_family or \
                    frozenset({axis_family, target_family}) in ADJACENT_PAIRS:
                adjacency_a2[j, t] = 1.0

    # ---- the two arms ---------------------------------------------------------------------
    matching = _match_on_legibility(table)
    index_of = {str(r["axis"]): int(r["axis_index"]) for _, r in table.iterrows()}
    certified_rows = [index_of[p["certified_axis"]] for p in matching["pairs"]]
    control_rows = [index_of[p["control_axis"]] for p in matching["pairs"]]

    def _arm(axis_rows, adjacency, label):
        sub_readout = readout[axis_rows]
        sub_adjacency = adjacency[axis_rows]
        finite = np.isfinite(sub_adjacency).all(axis=0) & np.isfinite(sub_readout).all(axis=0)
        sub_readout, sub_adjacency = sub_readout[:, finite], sub_adjacency[:, finite]
        if sub_readout.shape[1] < 3 or np.nanstd(sub_adjacency) < 1e-12:
            return {"arm": label, "status": "degenerate", "n_targets": int(finite.sum())}
        null = _permutation_null(sub_readout, sub_adjacency, n_draws=args.n_permutations,
                                 seed=args.seed)
        target_only = _spearman(
            (sub_readout - sub_readout.mean(axis=0, keepdims=True)).ravel(),
            (sub_adjacency - sub_adjacency.mean(axis=0, keepdims=True)).ravel())
        loo = {}
        for position, axis in enumerate(axis_rows):
            keep_rows = [i for i in range(len(axis_rows)) if i != position]
            loo[axis_name.get(axis, str(axis))] = _association(sub_readout[keep_rows],
                                                               sub_adjacency[keep_rows])
        return {"arm": label, "status": "scored", "n_axes": len(axis_rows),
                "n_targets": int(finite.sum()),
                "target_names_used": [str(n) for n in np.asarray(target_names)[finite]],
                **null, "target_centred_only_spearman": target_only,
                "leave_one_axis_out": loo,
                "leave_one_axis_out_min": float(min(loo.values())) if loo else float("nan"),
                "leave_one_axis_out_max": float(max(loo.values())) if loo else float("nan")}

    results = {}
    for label, adjacency in (("A1", adjacency_a1), ("A2", adjacency_a2)):
        certified = _arm(certified_rows, adjacency, f"certified_{label}")
        control = _arm(control_rows, adjacency, f"matched_uncertified_{label}")
        difference = {"status": "not_scored"}
        if certified.get("status") == "scored" and control.get("status") == "scored":
            rng = np.random.default_rng(args.seed)
            draws = []
            for _ in range(args.n_boot):
                pick = rng.integers(0, len(certified_rows), size=len(certified_rows))
                a_rows = [certified_rows[i] for i in pick]
                b_rows = [control_rows[i] for i in pick]
                finite = np.isfinite(adjacency[a_rows]).all(axis=0) \
                    & np.isfinite(adjacency[b_rows]).all(axis=0) \
                    & np.isfinite(readout[a_rows]).all(axis=0) \
                    & np.isfinite(readout[b_rows]).all(axis=0)
                if finite.sum() < 3:
                    continue
                draws.append(_association(readout[a_rows][:, finite],
                                          adjacency[a_rows][:, finite])
                             - _association(readout[b_rows][:, finite],
                                            adjacency[b_rows][:, finite]))
            draws = np.asarray(draws, dtype=float)
            draws = draws[np.isfinite(draws)]
            if draws.size:
                difference = {
                    "status": "scored", "n_draws": int(draws.size),
                    "point": float(certified["observed"] - control["observed"]),
                    "ci95_low": float(np.percentile(draws, 2.5)),
                    "ci95_high": float(np.percentile(draws, 97.5)),
                    "ci_excludes_zero": bool(np.percentile(draws, 2.5) > 0
                                             or np.percentile(draws, 97.5) < 0)}
        bridge = bool(certified.get("status") == "scored"
                      and certified["observed"] >= 0.15
                      and certified["permutation_p"] < 0.05
                      and difference.get("status") == "scored"
                      and difference["point"] > 0 and difference["ci_excludes_zero"])
        results[label] = {"certified": certified, "matched_uncertified": control,
                          "difference": difference, "bridge_works": bridge}

    # per-family breakdown of the certified arm, A2 (distrust check §2.8)
    per_family = {}
    for position, axis in enumerate(certified_rows):
        family = AXIS_FAMILY.get(axis_name.get(axis, ""), "?")
        adjacent = adjacency_a2[axis] > 0
        if adjacent.any() and (~adjacent).any():
            per_family.setdefault(family, []).append(
                float(np.nanmean(readout[axis][adjacent]) - np.nanmean(readout[axis][~adjacent])))
    per_family = {k: {"n_axes": len(v), "mean_adjacent_minus_non_adjacent": float(np.mean(v))}
                  for k, v in per_family.items()}

    n_significant = int(np.sum(np.asarray(adjacency_a1_lit)[np.isfinite(adjacency_a1_lit)] < 0.05))
    out = {
        "artifact": str(Path(args.artifact).name), "state": args.state,
        "partition": args.partition, "n_patients": int(len(rows)),
        "n_train": int(len(train)), "n_test": int(len(test)),
        "n_components": N_COMPONENTS, "caliper": CALIPER,
        "n_certified_axes": int((table["causal_name_certified"] == True).sum()),  # noqa: E712
        "targets": [{"name": str(n), "group": str(g), "family_a2": family_of_target[str(n)],
                     "has_gene_set": bool(str(n) in gene_sets)}
                    for n, g in zip(target_names, target_groups)],
        "targets_without_a_gene_set": missing_sets,
        "targets_unmapped_by_a2": unmapped,
        "matching": matching,
        "self_readout_median": float(np.nanmedian(self_readout)),
        "self_readout_vs_published_legibility_spearman": _spearman(
            self_readout[[int(r["axis_index"]) for _, r in table.iterrows()]],
            table[LEGIBILITY].to_numpy(float)),
        "a1_lit": {"n_pairs_scored": int(np.isfinite(adjacency_a1_lit).sum()),
                   "n_pairs_significant_p_lt_0_05": n_significant,
                   "underpowered": bool(n_significant < 10)},
        "results": results,
        "per_axis_family_a2": per_family,
        "readout_matrix": readout.tolist(),
        "adjacency_a1": adjacency_a1.tolist(),
        "adjacency_a2": adjacency_a2.tolist(),
        "axis_names": [axis_name.get(j, str(j)) for j in range(n_axes)],
        "perturbation_provenance": primary["provenance"],
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("matching", "a1_lit", "results", "per_axis_family_a2",
                                          "targets_without_a_gene_set")}, indent=2, default=float),
          flush=True)
    print(f"[done] -> {args.output}", flush=True)


if __name__ == "__main__":
    main()
