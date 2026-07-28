"""Evaluate closed-RAG Qwen cards against hidden held-out RNA targets only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .evidence_registry import EvidenceRegistry, default_evidence_registry
from .qwen_cards import evaluate_qwen_cards


def _records(path: str | Path) -> list[dict]:
    source = Path(path)
    content = source.read_text(encoding="utf-8").strip()
    if not content:
        raise ValueError(f"empty card file: {source}")
    if source.suffix.lower() == ".jsonl":
        values = [json.loads(line) for line in content.splitlines() if line.strip()]
    else:
        values = json.loads(content)
    if not isinstance(values, list):
        raise ValueError("card input must be JSONL or a JSON list")
    return values


def _table(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    return pd.read_parquet(source) if source.suffix.lower() in {".parquet", ".pq"} else pd.read_csv(source)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qwen-cards", required=True)
    parser.add_argument("--source-cards", required=True)
    parser.add_argument("--targets", required=True, help="RNA target table; never read by Qwen generation")
    parser.add_argument("--metadata", required=True, help="canonical patient metadata containing split")
    parser.add_argument("--output", required=True)
    parser.add_argument("--registry", default="")
    parser.add_argument("--k", type=int, default=3)
    args = parser.parse_args()
    metadata = _table(args.metadata)
    if not {"patient_id", "split"}.issubset(metadata):
        raise ValueError("metadata needs patient_id and split")
    heldout = set(metadata.loc[metadata.split.astype(str).eq("test"), "patient_id"].astype(str))
    targets = _table(args.targets)
    if "patient_id" not in targets:
        raise ValueError("targets need patient_id")
    targets["patient_id"] = targets.patient_id.astype(str)
    programmes = [column for column in targets.columns if column != "patient_id"]
    target_test = targets.loc[targets.patient_id.isin(heldout)].set_index("patient_id")
    source = [card for card in _records(args.source_cards) if str(card["patient_id"]) in heldout]
    qwen = [card for card in _records(args.qwen_cards) if str(card["patient_id"]) in heldout]
    ids = sorted(heldout & set(target_test.index.astype(str)))
    registry = EvidenceRegistry.from_json(args.registry) if args.registry else default_evidence_registry()
    metrics = evaluate_qwen_cards(
        qwen,
        source,
        ids,
        target_test.loc[ids, programmes].to_numpy(float),
        programmes,
        k=args.k,
        registry=registry,
    )
    payload = {"split": "test", "n_target_patients": len(ids), "registry_version": registry.version,
               "registry_digest": registry.digest, "metrics": metrics}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
