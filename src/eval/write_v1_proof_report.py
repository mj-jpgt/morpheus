"""Write a compact v1 proof report from emitted artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from morpheus.src.utils.config import load_config


def write_report(config_path: str = "morpheus/configs/v1.json") -> Path:
    cfg = load_config(config_path)
    out_dir = cfg.path("outputs_dir") / "v1_final"
    out_dir.mkdir(parents=True, exist_ok=True)
    audit_path = cfg.path("outputs_dir") / "v1_data_audit" / "data_inventory.json"
    align_dir = cfg.path("outputs_dir") / "v1_alignment"
    mp_dir = cfg.path("outputs_dir") / "v1_molecular_prompting"
    lines = ["# MORPHEUS V1 Proof Report", ""]
    lines.append("## Claim Level")
    lines.append("This report is predictive/retrieval-based proof only. It is not a clinical recommendation or validated biomarker report.")
    lines.append("")
    bulkformer_store = cfg.path("rna_bulkformer_embeddings")
    lines.append("## Foundation Model Status")
    if bulkformer_store.exists():
        lines.append(f"- RNA encoder: BulkFormer frozen embeddings from `{bulkformer_store}`.")
        extraction_manifest = bulkformer_store.with_name("tcga_bulkformer_extraction_manifest.json")
        if extraction_manifest.exists():
            manifest = json.load(open(extraction_manifest, encoding="utf-8"))
            lines.append(f"- BulkFormer extraction: n={manifest.get('n_samples')}, dim={manifest.get('embedding_dim')}, device={manifest.get('device')}.")
    else:
        lines.append("- RNA encoder: BulkFormer embeddings are not present; this is a fallback/smoke report, not the official v1 proof.")
    wsi_standard = cfg.path("wsi_standard_dir") / "tcga_ut_hoptimus0_manifest.json"
    if wsi_standard.exists():
        lines.append(f"- WSI encoder: H-Optimus-0 standardized feature store from `{wsi_standard}`.")
    else:
        lines.append("- WSI encoder: using legacy Morpheus WSI feature store if present.")
    lines.append("")
    if audit_path.exists():
        audit = json.load(open(audit_path, encoding="utf-8"))
        lines.append("## Data Readiness")
        for name, record in audit.get("modalities", {}).items():
            lines.append(f"- `{name}`: {record.get('counts', {})}")
        if audit.get("bulkformer_status", {}).get("official_v1_blocked_until_configured") and not bulkformer_store.exists():
            lines.append("- Blocking: BulkFormer repo ID is not configured; official v1 RNA encoder remains pending.")
        lines.append("")
    lines.append("## Alignment Metrics")
    if align_dir.exists():
        for metric_path in sorted(align_dir.glob("*_metrics.json")):
            payload = json.load(open(metric_path, encoding="utf-8"))
            metrics = payload.get("metrics", {})
            lines.append(
                f"- `{payload.get('method')}` ({payload.get('rna_source', 'unknown RNA source')}): "
                f"R@1={_fmt(metrics.get('recall_at_1'))}, R@5={_fmt(metrics.get('recall_at_5'))}, "
                f"R@10={_fmt(metrics.get('recall_at_10'))}, MRR={_fmt(metrics.get('mrr'))}, "
                f"sameCancer@10={_fmt(metrics.get('same_cancer_in_top10'))}, n_eval={payload.get('n_eval_pairs')}"
            )
    else:
        lines.append("- No alignment metrics found.")
    lines.append("")
    lines.append("## Molecular Prompting")
    mp_paths = sorted(mp_dir.glob("*/global_metrics.json")) if mp_dir.exists() else []
    if not mp_paths and (mp_dir / "global_metrics.json").exists():
        mp_paths = [mp_dir / "global_metrics.json"]
    if mp_paths:
        for path in mp_paths:
            mp = json.load(open(path, encoding="utf-8"))
            label = path.parent.name if path.parent != mp_dir else "default"
            lines.append(
                f"- `{label}` ({mp.get('target_source')}): mean R2={_fmt(mp.get('mean_r2'))}, "
                f"mean Pearson={_fmt(mp.get('mean_pearson'))}, mean Spearman={_fmt(mp.get('mean_spearman'))}, "
                f"gene sets={mp.get('n_gene_sets')}"
            )
    else:
        lines.append("- Not run.")
    lines.append("")
    lines.append("## Next Gate")
    lines.append("Proceed to v2 Query Former / multimodal expansion after preserving these frozen-encoder baselines, adding selected non-WSI omics tables, and deciding whether to keep CCA or CLIP as the v1 reference alignment.")
    report_path = out_dir / "v1_proof_report.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def _fmt(value) -> str:
    if value is None:
        return "NA"
    try:
        return f"{float(value):.4f}"
    except (TypeError, ValueError):
        return str(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="morpheus/configs/v1.json")
    args = parser.parse_args()
    print(write_report(args.config))


if __name__ == "__main__":
    main()
