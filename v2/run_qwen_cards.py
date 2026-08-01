"""CLI for the optional closed-RAG Qwen card-rendering stage.

Usage intentionally requires ``--enable-qwen`` before loading a model.  The
default emits validated abstentions, which keeps an unattended run safe when
weights are unavailable or the upstream quantitative gate has not passed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .card_quality_gate import load_quality_gate
from .evidence_registry import EvidenceRegistry, default_evidence_registry
from .qwen_cards import QwenGenerationConfig, TransformersQwenBackend, render_qwen_cards, write_qwen_jsonl


def _load_cards(path: Path) -> list[dict]:
    text = path.read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError("input card file is empty")
    if path.suffix.lower() == ".jsonl":
        values = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        values = json.loads(text)
    if not isinstance(values, list):
        raise ValueError("input must be a JSON array or JSONL sequence of cards")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description="Closed-RAG, schema-validated optional Qwen card rendering")
    parser.add_argument("--cards", required=True, type=Path, help="Deterministic hypothesis-card JSON or JSONL")
    parser.add_argument("--output", required=True, type=Path, help="Qwen-card JSONL output")
    parser.add_argument("--registry", type=Path, help="Versioned local evidence registry JSON; default is bundled ontology")
    parser.add_argument("--enable-qwen", action="store_true", help="Explicitly load local Qwen weights; default is safe abstention")
    parser.add_argument("--model-id", default="Qwen/Qwen2.5-7B-Instruct")
    parser.add_argument("--allow-download", action="store_true", help="Permit Hugging Face model download; disabled by default")
    parser.add_argument("--quality-gate", type=Path, help="required passed quantitative card-gate JSON when enabling Qwen")
    parser.add_argument("--max-new-tokens", type=int, default=384)
    args = parser.parse_args()
    if args.enable_qwen:
        if args.quality_gate is None:
            raise ValueError("--enable-qwen requires --quality-gate")
        gate = load_quality_gate(args.quality_gate)
        if not gate["passed"]:
            raise ValueError("Qwen rendering is disabled because the quantitative card gate did not pass")
    else:
        gate = {"passed": False, "reason": "qwen_not_enabled"}
    cards = _load_cards(args.cards)
    registry = EvidenceRegistry.from_json(args.registry) if args.registry else default_evidence_registry()
    config = QwenGenerationConfig(model_id=args.model_id, max_new_tokens=args.max_new_tokens, local_files_only=not args.allow_download)
    backend = TransformersQwenBackend(config) if args.enable_qwen else None
    results = render_qwen_cards(cards, registry, backend, config=config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_qwen_jsonl(results, args.output)
    manifest = {
        "cards": str(args.cards), "output": str(args.output), "registry_version": registry.version,
        "registry_digest": registry.digest, "qwen_enabled": bool(args.enable_qwen), "model_id": args.model_id if args.enable_qwen else None,
        "n_cards": len(results), "quality_gate": gate, "quality_gate_path": str(args.quality_gate) if args.quality_gate else None,
    }
    args.output.with_suffix(args.output.suffix + ".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
