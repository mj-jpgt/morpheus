"""Write diagnostics comparing MORPHEUS V2 representations with baseline summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from morpheus.src.eval.retrieval_metrics import l2_normalize, paired_retrieval_metrics


DEFAULT_EMBEDDINGS = {
    "v2_patient_stratified": Path("morpheus/outputs/v2_bio_query_former_patient_stratified_full/bio_query_former_embeddings.npz"),
    "v2_patient_hybrid_fn": Path("morpheus/outputs/v2_bio_query_former_patient_stratified_hybrid_fn_full/bio_query_former_embeddings.npz"),
    "v2_patch_stratified": Path("morpheus/outputs/v2_bio_query_former_patch_stratified_full/bio_query_former_embeddings.npz"),
}


def _mean_pairwise_same_cancer_top10(wsi: np.ndarray, rna: np.ndarray, cancers: np.ndarray) -> float:
    metrics = paired_retrieval_metrics(wsi, rna, (10,), cancers.astype(str).tolist(), cancers.astype(str).tolist())
    return float(metrics.get("same_cancer_in_top10", 0.0))


def _mean_diag_cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = l2_normalize(a)
    b = l2_normalize(b)
    return float(np.mean(np.sum(a * b, axis=1)))


def _cancer_leakage(x: np.ndarray, cancers: np.ndarray, split: np.ndarray, max_train: int = 1500, max_test: int = 1500) -> float | None:
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import balanced_accuracy_score
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import LabelEncoder, StandardScaler

    train = split == "train"
    test = split == "test"
    if train.sum() < 10 or test.sum() < 10 or len(np.unique(cancers[train])) < 2:
        return None
    rng = np.random.default_rng(42)
    train_idx = np.where(train)[0]
    test_idx = np.where(test)[0]
    if len(train_idx) > max_train:
        train_idx = np.sort(rng.choice(train_idx, size=max_train, replace=False))
    if len(test_idx) > max_test:
        test_idx = np.sort(rng.choice(test_idx, size=max_test, replace=False))
    train = np.zeros(len(split), dtype=bool)
    test = np.zeros(len(split), dtype=bool)
    train[train_idx] = True
    test[test_idx] = True
    encoder = LabelEncoder()
    y_train = encoder.fit_transform(cancers[train])
    known = np.isin(cancers[test], encoder.classes_)
    if known.sum() < 10:
        return None
    model = make_pipeline(StandardScaler(), OneVsRestClassifier(LogisticRegression(max_iter=500, class_weight="balanced", solver="liblinear")))
    model.fit(x[train], y_train)
    pred = model.predict(x[test][known])
    return float(balanced_accuracy_score(encoder.transform(cancers[test][known]), pred))


def _embedding_diagnostics(name: str, path: Path) -> list[dict]:
    data = np.load(path, allow_pickle=True)
    split = data["split"].astype(str)
    cancers = data["cancers"].astype(str)
    rows = []
    for feature in ("wsi_identity", "wsi_biology", "full_biology"):
        if feature not in data.files:
            continue
        x = data[feature].astype(np.float32)
        row = {
            "run": name,
            "feature": feature,
            "n": int(len(x)),
            "test_cancer_balanced_accuracy": _cancer_leakage(x, cancers, split),
        }
        if feature != "full_biology" and "rna_identity" in data.files:
            row["same_cancer_at_10"] = _mean_pairwise_same_cancer_top10(data["wsi_identity"], data["rna_identity"], cancers)
        if "wsi_identity" in data.files and feature != "wsi_identity":
            row["identity_feature_diag_cosine"] = _mean_diag_cosine(data["wsi_identity"], x)
        rows.append(row)
    if "wsi_program_scores" in data.files and data["wsi_program_scores"].shape[1] > 0:
        rows.append(
            {
                "run": name,
                "feature": "wsi_program_scores",
                "n": int(len(data["wsi_program_scores"])),
                "test_cancer_balanced_accuracy": _cancer_leakage(data["wsi_program_scores"].astype(np.float32), cancers, split),
            }
        )
    return rows


def _read_optional_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return ""
    values = frame.fillna("NA").astype(str)
    header = "| " + " | ".join(values.columns) + " |"
    sep = "| " + " | ".join("---" for _ in values.columns) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in values.to_numpy()]
    return "\n".join([header, sep, *rows])


def write_v2_diagnostics_report(output_dir: str | Path = "morpheus/outputs/v2_diagnostics") -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, path in DEFAULT_EMBEDDINGS.items():
        if path.exists():
            rows.extend(_embedding_diagnostics(name, path))
    diag = pd.DataFrame(rows)
    diag.to_csv(out / "representation_diagnostics.csv", index=False)

    v2_molecular = _read_optional_csv(Path("morpheus/outputs/v2_molecular_benchmarks/v2_molecular_benchmark_summary.csv"))
    fair_molecular = _read_optional_csv(Path("morpheus/outputs/fair_molecular_benchmarks/fair_molecular_benchmark_summary.csv"))
    if not v2_molecular.empty:
        v2_molecular.to_csv(out / "v2_molecular_summary_snapshot.csv", index=False)
    if not fair_molecular.empty:
        fair_molecular.to_csv(out / "fair_molecular_summary_snapshot.csv", index=False)

    lines = [
        "# V2 Diagnostics Report",
        "",
        "Diagnostics are computed on fixed saved splits and do not refit or alter train/val/test assignments.",
        "",
        "## Representation Diagnostics",
        "",
        _markdown_table(diag) if not diag.empty else "No V2 embedding diagnostics were available.",
        "",
        "## Molecular Transfer Summary",
        "",
    ]
    if not v2_molecular.empty:
        cols = [c for c in ["run", "feature", "test_soft_knn_pearson", "test_common_head_pearson", "test_n"] if c in v2_molecular.columns]
        lines.extend(["### V2", "", _markdown_table(v2_molecular[cols]), ""])
    if not fair_molecular.empty:
        cols = [c for c in ["split", "method", "test_soft_knn_pearson", "test_common_head_pearson", "test_n"] if c in fair_molecular.columns]
        lines.extend(["### Fair Baselines", "", _markdown_table(fair_molecular[cols]), ""])
    lines.extend(
        [
            "## Interpretation Guardrails",
            "",
            "- High cancer balanced accuracy means the feature still carries lineage/cancer identity and should not be claimed as cancer-invariant biology.",
            "- `full_biology` includes fused multimodal context and should not be compared as WSI-only molecular prompting.",
            "- Patch-mode rows using local 2048-dimensional patch bags are not H-Optimus-0 patch-token evidence until provenance is confirmed or re-extracted.",
        ]
    )
    report = out / "v2_diagnostics_report.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "v2_diagnostics_manifest.json").write_text(json.dumps({"embedding_jobs": {k: str(v) for k, v in DEFAULT_EMBEDDINGS.items()}}, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="morpheus/outputs/v2_diagnostics")
    args = parser.parse_args()
    print(write_v2_diagnostics_report(args.output_dir))


if __name__ == "__main__":
    main()
