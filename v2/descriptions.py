"""Read and audit the public cancer-description resource used by V2."""
from __future__ import annotations
from pathlib import Path
from .text_prototypes import CancerDescription, validate_descriptions


def load_tcga_descriptions(path: str | Path) -> list[CancerDescription]:
    try:
        import yaml
    except ImportError as exc:
        raise RuntimeError("PyYAML is required to read cancer descriptions") from exc
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    source = raw["source_url"]
    values = []
    for cancer, item in raw["cancers"].items():
        text = ". ".join((item["name"], "Tissue of origin: " + item["tissue"], "Histology: " + item["histology"], "Canonical features: " + item["features"], "Context: " + item["context"]))
        values.append(CancerDescription(cancer=str(cancer), text=text, source=source, public_knowledge_only=True))
    validate_descriptions(values, raw["cancers"].keys())
    return values
