"""P2: competing label-free representation-quality metrics, evaluated AS SELECTION RULES.

VENDORED 2026-08-04 from `~/e0_run/p2_competing_metrics.py` on `ubuntu@150.136.45.194`, the
script that produced `paper/P2_RANK_DRAFT.md` §4.2, §4.4 and §4.6. It is in the repository
because §6.2 lists "the competing-metric and variance-decomposition scripts, vendored into the
repository" as a pre-submission requirement: the paper's own standard is that every number
traces to a file in the repository, and until this commit these four scripts existed only on
the box. The vendored copy is byte-for-byte the original except for this paragraph and for
`participation_ratio_rownorm`, added below and used by no published number.

Provenance of the reproduction: `NOTEBOOK_ENTRIES/p2_vendored_and_reproduced_20260804T0255Z.md`.

Computes, on the frozen D2 and D1 representation artifacts:
  * canonical effective rank (Roy-Vetterli, this repository's single source of truth)
  * RankMe (Garrido, Balestriero, Najman & LeCun, ICML 2023, arXiv:2210.02885)
  * alpha-ReQ  [filled in only once the estimator is verified from source]
  * LiDAR (Thilak et al., arXiv:2312.04000)  [ditto]

and the ground-truth molecular channel used by d2_readout.py, so that each metric can be
scored on the only question that matters for the paper: DOES IT PICK THE ARM THAT ACTUALLY
CARRIES MORE MOLECULAR INFORMATION?

Ground truth is the top canonical correlation between the residualised representation and
the residualised 40 untrained targets, computed with EXACTLY the residualisation, seed and
component budget of d2_readout.py so the numbers are comparable to D2_PER_ARTIFACT_READOUT.json.

METHODOLOGICAL NOTE, PRE-DECLARED. Every metric is computed twice: once on the raw
representation and once on the confound-residualised representation. Effective rank in this
project is reported on the residualised representation, but RankMe as published is computed
on the raw embedding matrix. Reporting only one would let a referee argue that the comparison
was rigged by the preprocessing, so both are reported and both are used to score the
selection rule. If the two disagree, that disagreement is itself a result.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.v2.calibra.residualise import (confound_design, cross_fitted_residuals,
                                             pooled_tissue_source_site)
from morpheus.v2.calibra.spectral import effective_rank, top_canonical_correlation

GROUPS = {
    "unrestricted": None,
    "untrained40": ["heldout_pathway", "immune_tme", "tumour_state"],
    "random_control": ["random_control"],
    "hallmark_in_training": ["hallmark_in_training"],
}


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def rankme(z, eps: float = 1e-7) -> float:
    """RankMe, exactly as published.

    Garrido, Balestriero, Najman & LeCun, "RankMe: Assessing the Downstream Performance of
    Pretrained Self-Supervised Representations by Their Rank", ICML 2023, arXiv:2210.02885.

    The published definition is Roy-Vetterli effective rank with an epsilon inside the
    normalisation:  p_k = sigma_k(Z) / ||sigma(Z)||_1  +  eps , then exp(-sum p_k log p_k).

    Two deliberate differences from this repository's `effective_rank`, both of which follow
    the published formula rather than our convention:
      (1) NO column centring. RankMe is defined on the embedding matrix Z as it comes out of
          the network. Our effective_rank column-centres first.
      (2) NO zero-singular-value pruning; the epsilon is what keeps log finite.
    """
    z = np.asarray(z, dtype=np.float64)
    singular = np.linalg.svd(z, compute_uv=False)
    p = singular / np.abs(singular).sum() + eps
    return float(np.exp(-(p * np.log(p)).sum()))


def spectrum(x, centre: bool = True) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if centre:
        x = x - x.mean(axis=0, keepdims=True)
    return np.linalg.svd(x, compute_uv=False)


def participation_ratio(x) -> float:
    """PR = (sum lambda)^2 / sum(lambda^2) on the eigenvalues of the covariance.

    Included as a second, independent linear-dimensionality proxy so that the paper's claim
    is not resting on one particular functional form of "effective rank". Not attributed to
    any specific paper here; it is a standard second-moment dimensionality statistic.

    NOT `RANK_VARIANTS["R2"]`, despite both being called a participation ratio. This one is
    the order-2 Hill number of the EIGENVALUE distribution, `(sum s^2)^2 / sum s^4`; R2 is
    the order-2 Hill number of the SINGULAR-VALUE distribution, `(sum s)^2 / sum s^2`, which
    is what `d1_audit.py` computed and what `calibra/spectral.py` now registers under that
    name. The two are different functions of the same spectrum and R2 is the larger of the
    two for every matrix. Which of them a "participation ratio" row was computed with is a
    thing the paper has to state; see `p2_rank_variants.py`, which reports both.
    """
    lam = spectrum(x) ** 2
    return float(lam.sum() ** 2 / (lam ** 2).sum())


def participation_ratio_rownorm(x) -> float:
    """`participation_ratio` on L2-normalised rows -- the eigenvalue analogue of R3.

    ADDED AT VENDORING TIME and used by no published number. It exists only so that
    `p2_rank_variants.py` can report the statistic its pre-vendoring version actually
    computed for its "R3" column beside the canonical `RANK_VARIANTS["R3"]`, rather than
    asserting that the two agree.
    """
    value = np.asarray(x, dtype=np.float64)
    norms = np.linalg.norm(value, axis=1, keepdims=True)
    norms[norms < 1e-12] = 1.0
    return participation_ratio(value / norms)


def stable_rank(x) -> float:
    """||X||_F^2 / ||X||_2^2 -- the standard stable rank. Same purpose as above."""
    s = spectrum(x)
    return float((s ** 2).sum() / s[0] ** 2)


# --- alpha-ReQ ------------------------------------------------------------
def alpha_req(x, trange_lo: int = 3, trange_hi: int = 100, centre: bool = True):
    """alpha-ReQ eigenspectrum-decay exponent, following the AUTHORS' RELEASED CODE.

    Agrawal, Mondal, Ghosh & Richards, "alpha-ReQ: Assessing Representation Quality in
    Self-Supervised Learning by measuring eigenspectrum decay", NeurIPS 2022, vol. 35,
    pp. 17626-17638, doi:10.52202/068431-1281.

    The PAPER says only "we compute the full set of numerical eigenvalues, and estimate alpha
    by fitting a linear model on the eigenspectrum in log-log scale" and gives NO index range.
    The index range and the fit weighting exist only in the authors' `fastssl` repository
    (`fastssl/utils/powerlaw.py`, `stringer_get_powerlaw`, called from `scripts/train_model.py`
    as `trange=np.arange(3, 100)`). This reimplementation follows the CODE, because the paper
    text is not sufficient to reproduce a number:

      * centre the features, PCA, take `explained_variance_ratio_` as the spectrum;
      * weighted least squares of log(lambda) on -log(rank) over eigenvalue ranks 4..100
        (0-based indices 3..99), with weights w = 1/rank.

    The discrepancy (paper says "full set", code fits ranks ~4-100) is real and is recorded in
    the notebook entry; the index range must NOT be attributed to the paper text.
    """
    x = np.asarray(x, dtype=np.float64)
    if centre:
        x = x - x.mean(axis=0, keepdims=True)
    s_ = np.linalg.svd(x, compute_uv=False)
    lam = s_ ** 2
    lam = lam / lam.sum()                      # explained_variance_ratio_ equivalent
    trange = np.arange(trange_lo, min(trange_hi, len(lam)))
    if len(trange) < 5:
        return {"alpha": float("nan"), "n_points": 0}
    logss = np.log(np.abs(lam))
    y = logss[trange][:, np.newaxis]
    tr = trange + 1                            # authors increment AFTER indexing
    nt = tr.size
    xd = np.concatenate((-np.log(tr)[:, np.newaxis], np.ones((nt, 1))), axis=1)
    w = 1.0 / tr.astype(np.float64)[:, np.newaxis]
    b = np.linalg.solve(xd.T @ (xd * w), (w * xd).T @ y).flatten()
    pred = (xd @ b[:, None]).ravel()
    resid = y.ravel() - pred
    r2 = 1.0 - resid.var() / y.ravel().var() if y.ravel().var() > 0 else float("nan")
    return {"alpha": float(b[0]), "n_points": int(nt), "fit_r2": float(r2),
            "trange": [int(trange_lo), int(trange[-1] + 1)],
            "source": "authors fastssl/utils/powerlaw.py stringer_get_powerlaw, trange=arange(3,100)"}


# --- LiDAR ----------------------------------------------------------------
def lidar(views, delta: float = 1e-4, eps: float = 1e-7):
    """LiDAR: linear-discriminant rank of a joint-embedding representation.

    Thilak, Huang, Saremi, Dinh, Goh, Nakkiran, Susskind & Littwin, arXiv:2312.04000.

    `views` is a list of (n, d) arrays: view v of the SAME n objects. Each object is one
    class of the surrogate discriminative task; the views are the within-class samples.

    ADAPTATION, DECLARED -- but explicitly licensed by the paper. LiDAR footnote 4 defines the
    within-class samples as "Data augmentations, OR OTHERWISE DATA POINTS WHICH ARE TREATED AS
    POSITIVE SAMPLES". Our joint-embedding objective treats a patient's WSI-derived and
    RNA-derived embeddings as exactly that positive pair, so using the two modalities as the two
    views is sanctioned by the definition. What is NOT in the authors' tested regime is q:
    they use q = 50 (recommended) or q = 10, sweep q only over 10..100, and never discuss q = 2.
    Their n > p recommendation IS met (n = 2766 patients > p = 256 dims).

    LiDAR's surrogate task otherwise uses q augmentations of each clean sample.
    Our frozen artifacts do not store augmentations -- they store the model's per-patient
    embedding under each MODALITY (WSI-derived and RNA-derived) of the joint-embedding
    model. We therefore use the modalities as the views, with q = 2, one patient = one
    class. This is a faithful use of the surrogate-task construction (invariance directions
    are exactly what the joint-embedding objective is trained to produce) but it is NOT the
    augmentation-based instantiation the authors benchmarked, and the paper must say so.
    q = 2 is the minimum that makes the within-class covariance estimable at all.

    Sigma_w = mean over classes of the within-class scatter, + delta * I
    Sigma_b = covariance of the class means
    LiDAR matrix = Sigma_w^{-1/2} Sigma_b Sigma_w^{-1/2}; score = smooth rank of its
    eigenvalues (same entropy-exponential form as RankMe).
    """
    v = np.stack([np.asarray(a, dtype=np.float64) for a in views], axis=0)  # (q, n, d)
    q, n, d = v.shape
    class_mean = v.mean(axis=0)                       # (n, d)
    grand = class_mean.mean(axis=0, keepdims=True)    # (1, d)

    dev = (v - class_mean[None]).reshape(q * n, d)
    sigma_w = dev.T @ dev / (n * (q - 1))
    sigma_w = sigma_w + delta * np.eye(d)

    mb = class_mean - grand
    sigma_b = mb.T @ mb / (n - 1)

    w, u = np.linalg.eigh(sigma_w)
    w = np.clip(w, 1e-12, None)
    inv_sqrt = (u / np.sqrt(w)) @ u.T
    lidar_mat = inv_sqrt @ sigma_b @ inv_sqrt
    lam = np.linalg.eigvalsh(lidar_mat)
    lam = np.clip(lam, 0.0, None)
    if lam.sum() <= 0:
        return {"lidar": float("nan")}
    p = lam / lam.sum() + eps
    return {"lidar": float(np.exp(-(p * np.log(p)).sum())),
            "delta": delta, "q": int(q), "n_classes": int(n)}


# ---------------------------------------------------------------------------
def load_targets(path, groups):
    raw = np.load(path, allow_pickle=False)
    ids = np.asarray(raw["patient_ids"]).astype(str)
    names = np.asarray(raw["target_names"]).astype(str)
    values = np.asarray(raw["scores"], dtype=np.float64)
    available = np.asarray(raw["target_groups"]).astype(str)
    if groups:
        keep = np.isin(available, list(groups))
    else:
        keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
    return {p: i for i, p in enumerate(ids)}, values[:, keep], names[keep]


def score_artifact(path, targets_path, seed, alpha_range, lidar_delta, subsamples, sub_frac):
    raw = np.load(path, allow_pickle=False)
    test = raw["split"].astype(str) == "test"
    ids = raw["patient_ids"].astype(str)[test]
    cancers = raw["cancers"].astype(str)[test]
    x = raw["wsi_biology"].astype(np.float64)[test]
    rna = raw["rna_biology"].astype(np.float64)[test]

    tss, sites = pooled_tissue_source_site(ids, min_site_count=10)
    design = confound_design(pd.DataFrame({"cancer": cancers, "tss": tss}), ["cancer", "tss"])
    rx = cross_fitted_residuals(x, design, seed=seed)
    rrna = cross_fitted_residuals(rna, design, seed=seed)

    entry = {"path": str(Path(path).resolve()), "n_test": int(test.sum()),
             "metrics": {}, "points": {}}
    m = entry["metrics"]
    m["effective_rank_raw"] = effective_rank(x)
    m["effective_rank_residualised"] = effective_rank(rx)
    m["rankme_raw"] = rankme(x)
    m["rankme_residualised"] = rankme(rx)
    m["participation_ratio_raw"] = participation_ratio(x)
    m["participation_ratio_residualised"] = participation_ratio(rx)
    m["stable_rank_raw"] = stable_rank(x)
    m["stable_rank_residualised"] = stable_rank(rx)
    m["alpha_req_raw"] = alpha_req(x, *alpha_range)
    m["alpha_req_residualised"] = alpha_req(rx, *alpha_range)
    m["lidar_raw"] = lidar([x, rna], delta=lidar_delta)
    m["lidar_residualised"] = lidar([rx, rrna], delta=lidar_delta)

    # reproducibility floor: recompute the label-free metrics on patient subsamples so that
    # a between-arm gap can be compared against the metric's own sampling noise.
    if subsamples:
        rng = np.random.default_rng(seed)
        k = int(len(x) * sub_frac)
        acc = {"effective_rank_residualised": [], "rankme_raw": [],
               "participation_ratio_residualised": [], "lidar_residualised": []}
        for _ in range(subsamples):
            sel = rng.choice(len(x), size=k, replace=False)
            acc["effective_rank_residualised"].append(effective_rank(rx[sel]))
            acc["rankme_raw"].append(rankme(x[sel]))
            acc["participation_ratio_residualised"].append(participation_ratio(rx[sel]))
            acc["lidar_residualised"].append(lidar([rx[sel], rrna[sel]], delta=lidar_delta)["lidar"])
        entry["subsample"] = {kk: {"mean": float(np.mean(vv)), "sd": float(np.std(vv, ddof=1)),
                                   "n": len(vv), "frac": sub_frac}
                              for kk, vv in acc.items()}

    for gname, groups in GROUPS.items():
        index, values, _ = load_targets(targets_path, groups)
        rows = np.asarray([index[p] for p in ids])
        ry = cross_fitted_residuals(values[rows], design, seed=seed)
        entry["points"][gname] = {"n_targets": int(values.shape[1]),
                                  "top_cca": top_canonical_correlation(rx, ry, n_components=16)}
    return entry


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--artifacts", nargs="+", required=True, help="label=path")
    ap.add_argument("--targets", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alpha-index-range", nargs=2, type=int, default=[11, 50])
    ap.add_argument("--lidar-delta", type=float, default=1e-4)
    ap.add_argument("--subsamples", type=int, default=0)
    ap.add_argument("--subsample-fraction", type=float, default=0.8)
    a = ap.parse_args(argv)

    out = {"_config": {"seed": a.seed, "alpha_index_range": a.alpha_index_range,
                       "lidar_delta": a.lidar_delta, "subsamples": a.subsamples,
                       "subsample_fraction": a.subsample_fraction,
                       "targets": str(Path(a.targets).resolve())}}
    for item in a.artifacts:
        label, path = item.split("=", 1)
        entry = score_artifact(path, a.targets, a.seed, a.alpha_index_range, a.lidar_delta,
                               a.subsamples, a.subsample_fraction)
        out[label] = entry
        m = entry["metrics"]
        print(f"{label:24s} er_res={m['effective_rank_residualised']:8.3f} "
              f"rankme_raw={m['rankme_raw']:8.3f} "
              f"pr_res={m['participation_ratio_residualised']:8.3f} "
              f"alpha_res={m['alpha_req_residualised']['alpha']:6.3f} "
              f"lidar_res={m['lidar_residualised']['lidar']:8.3f} "
              f"cca40={entry['points']['untrained40']['top_cca']:.4f}", flush=True)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("wrote", a.output)
    return out


if __name__ == "__main__":
    main()
