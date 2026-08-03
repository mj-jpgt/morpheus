"""Is the programme_free biology head collapsing during REAL training?

Training metrics suggested it: arm F's per-feature std falls 0.034 -> 0.018 while
arm P's rises to 0.049 (isotropic for d=256 is 0.0625), and F's contrastive term
drifts from 6.15 toward its 8.37 chance value. But that metric is measured against
a LIVE queue that tracks the model, so it is confounded.

This measures the geometry directly on HELD-OUT test patients from the live
checkpoints: patient-to-patient cosine and effective rank of z_biology. Those are
unconfounded, and are what audit check A5 reports at the end.

usage: geom_probe.py <n_patients>
"""
import sys
import numpy as np, torch, torch.nn.functional as F
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.preflight import restrict_cohort_to_split
from morpheus.v2.runner import attach_v2_targets, UncappedHoptimusBatches, _truncate_batch
from morpheus.v2.model import TumorStateV2, V2ModelConfig

N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
D = "/home/ubuntu/e0_run/data/"
S = D + "paired_split_maximal.json"
torch.manual_seed(0)
data = load_bio_query_data(D + "v1_abs_hallmark.json", S, wsi_mode="hoptimus_patch")
data, _ = restrict_cohort_to_split(data, S)
fit = np.asarray(data.split).astype(str) != "test"
attach_v2_targets(data, fit)
data._v2_programme_head_dim = 256
test = np.where(np.asarray(data.split).astype(str) == "test")[0]
loader = UncappedHoptimusBatches(data, test, 8192, 0, shuffle=True)
batch = _truncate_batch(next(iter(loader)), N)
batch = {k: (v.cuda() if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
n = len(batch["indices"])
eye = torch.eye(n, dtype=torch.bool, device="cuda")
print(f"held-out test patients in probe: {n}\n")
print(f"{'arm':<20}{'epoch':>6}{'wsi-wsi cos':>13}{'rna-rna cos':>13}{'eff-rank':>10}{'feat-std':>10}{'rank@1e-3':>11}")
import os
ARMS = os.environ.get("ARMS", "p:programme_only:42,f:programme_free:42")
for spec in ARMS.split(","):
    arm, tag, sd = spec.split(":")
    ck = torch.load(f"/home/ubuntu/e0_run/d1_v1/d1_{arm}_seed{sd}/last.pt", map_location="cuda", weights_only=False)
    cfg = V2ModelConfig(**ck["manifest"]["model_config"])
    model = TumorStateV2(cfg, programme_dim=256).cuda()
    model.load_state_dict(ck["model"])
    model.eval()
    with torch.no_grad():
        w = model(batch, view="wsi")["z_biology"].float()
        r = model(batch, view="rna")["z_biology"].float()
    wn, rn = F.normalize(w, dim=-1), F.normalize(r, dim=-1)
    sv = torch.linalg.svdvals(wn - wn.mean(0))
    erank = float((sv.sum() ** 2) / (sv ** 2).sum())
    hard = int(torch.linalg.matrix_rank(wn.cpu(), tol=1e-3))
    print(f"{tag+chr(32)+chr(115)+sd:<20}{int(ck['epoch']):>6}{float((wn@wn.T)[~eye].mean()):>13.4f}"
          f"{float((rn@rn.T)[~eye].mean()):>13.4f}{erank:>10.2f}{float(wn.std(0).mean()):>10.4f}{hard:>11d}")
    del model, ck
    torch.cuda.empty_cache()
print("\nisotropic per-feature std for d=256 is 0.0625; eff-rank is the participation")
print("ratio of the centred singular values (max = n_patients).")
