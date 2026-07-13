"""Published and same-data baseline registry for TCGA-UT benchmark reports."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class BaselineSpec:
    name: str
    category: str
    input: str
    evaluation_role: str
    weights: str
    weights_url: str
    access: str
    notes: str


BASELINES: tuple[BaselineSpec, ...] = (
    BaselineSpec("H-Optimus-0 patch mean", "same encoder", "H-Optimus-0 tiles", "required", "H-Optimus-0", "https://huggingface.co/bioptimus/H-optimus-0", "gated", "Patch-weighted patient mean."),
    BaselineSpec("H-Optimus-0 slide-balanced mean", "same encoder", "H-Optimus-0 tiles", "required", "H-Optimus-0", "https://huggingface.co/bioptimus/H-optimus-0", "gated", "Mean each slide, then patient mean."),
    BaselineSpec("ABMIL on H-Optimus-0", "MIL", "frozen H-Optimus-0 tiles", "required", "train from protocol train", "", "train", "Attention MIL pooling with identical train/val/test protocol."),
    BaselineSpec("CLAM on H-Optimus-0", "MIL", "frozen H-Optimus-0 tiles", "required", "train from protocol train", "https://github.com/mahmoodlab/CLAM", "code public", "Cluster-constrained attention MIL."),
    BaselineSpec("Legacy patient-average BioQueryFormer", "ablation", "legacy patient means", "required", "local checkpoint", "", "local", "Historical comparator; cannot use token-level evidence."),
    BaselineSpec("Token-aware BioQueryFormer", "proposed", "H-Optimus-0 tokens", "primary", "train from protocol train", "", "train", "Slot aggregation and disentangled objectives."),
    BaselineSpec("Soft-kNN molecular prompting", "molecular", "aligned WSI/RNA", "required", "none", "", "none", "RNA reference set is train-only."),
    BaselineSpec("Ridge Hallmark head", "molecular", "WSI embeddings", "required", "none", "", "none", "StandardScaler plus RidgeCV fit train-only."),
    BaselineSpec("Phikon-v2", "foundation encoder", "tiles", "strong control", "Phikon-v2", "https://huggingface.co/owkin/phikon-v2", "public", "Use identical tile/pooling protocol."),
    BaselineSpec("UNI2-h", "foundation encoder", "tiles", "strong control", "UNI2-h", "https://huggingface.co/MahmoodLab/UNI2-h", "gated", "Academic/gated access; report revision and license."),
    BaselineSpec("CONCH", "vision-language encoder", "tiles", "retrieval control", "CONCH", "https://huggingface.co/MahmoodLab/CONCH", "gated", "Frozen image encoder with matched pooling."),
    BaselineSpec("Prov-GigaPath", "foundation encoder", "tiles", "strong control", "Prov-GigaPath", "https://huggingface.co/prov-gigapath/prov-gigapath", "gated", "Use only after input compatibility is recorded."),
    BaselineSpec("PathomicFusion/PORPOISE", "multimodal", "WSI plus omics", "optional task comparator", "train from protocol train", "https://github.com/mahmoodlab/PORPOISE", "code public", "Only compare on matching modalities and endpoints."),
    BaselineSpec("SurvPath", "multimodal", "WSI plus pathways", "optional task comparator", "train from protocol train", "https://github.com/mahmoodlab/SurvPath", "code public", "Survival-specific; not a retrieval baseline."),
)


def baseline_registry_frame() -> pd.DataFrame:
    return pd.DataFrame([asdict(spec) for spec in BASELINES])


def write_baseline_registry(output_dir: str | Path) -> tuple[Path, Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    table = baseline_registry_frame()
    csv_path = output / "tcga_ut_baseline_registry.csv"
    markdown_path = output / "tcga_ut_baseline_registry.md"
    table.to_csv(csv_path, index=False)
    markdown_path.write_text("# TCGA-UT Baseline Registry\n\n" + table.to_markdown(index=False) + "\n", encoding="utf-8")
    return csv_path, markdown_path
