"""Export frozen V2 representations in the shared patient-artifact schema."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import torch
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from .model import TumorStateV2, V2ModelConfig
from .runner import UncappedHoptimusBatches, attach_v2_targets


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data-config", required=True); p.add_argument("--split-file", required=True)
    p.add_argument("--checkpoint", required=True); p.add_argument("--output", required=True)
    p.add_argument("--token-budget", type=int, default=32768); p.add_argument("--hidden-dim", type=int, default=512)
    p.add_argument("--layers", type=int, default=4); p.add_argument("--heads", type=int, default=8); p.add_argument("--device", default="cuda")
    a = p.parse_args(); data = load_bio_query_data(a.data_config, a.split_file, wsi_mode="hoptimus_patch"); attach_v2_targets(data)
    cfg = V2ModelConfig(rna_dim=int(data.rna.shape[1]), hidden_dim=a.hidden_dim, layers=a.layers, heads=a.heads)
    model = TumorStateV2(cfg, programme_dim=int(data.hallmark.shape[1])).to(a.device)
    model.load_state_dict(torch.load(a.checkpoint, map_location=a.device, weights_only=False)["model"]); model.eval()
    n = len(data.patient_ids); arrays = {key: np.zeros((n, width), np.float32) for key, width in {"wsi_identity":256,"rna_identity":256,"wsi_biology":256,"full_identity":256,"full_biology":256,"full_patient":256,"uncertainty":cfg.hidden_dim}.items()}
    loader = UncappedHoptimusBatches(data, np.arange(n), a.token_budget, 917)
    with torch.no_grad():
        for batch in loader:
            indices = batch.pop("indices").numpy()
            batch = {k: v.to(a.device, non_blocking=True) for k,v in batch.items()}
            wsi, rna, full = model(batch,"wsi"), model(batch,"rna"), model(batch,"full")
            values = {"wsi_identity":wsi["z_identity"],"rna_identity":rna["z_identity"],"wsi_biology":wsi["z_biology"],"full_identity":full["z_identity"],"full_biology":full["z_biology"],"full_patient":full["z_patient"],"uncertainty":full["z_uncertainty"]}
            for key,value in values.items(): arrays[key][indices] = value.cpu().float().numpy()
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.output, patient_ids=np.asarray(data.patient_ids), split=np.asarray(data.split), cancers=np.asarray(data.cancers), **arrays)


if __name__ == "__main__": main()
