"""Does `decorrelation` cause programme_free's collapse at TRAINING scale?

The existing causal evidence is all from the G2.6 memorisation check: one fixed
16-patient batch, a frozen queue, no streaming. The observed D1 collapse is on
real streaming batches with a live 4096-key queue. This closes that gap.

Two throwaway models, identical in every respect except `decorrelation_after_warmup`
(0.04, as D1 runs, vs 0.0). Real train batches, live queue, bf16 autocast and
grad-clip 1.0 exactly as `V2Trainer.train_epoch` does. Centred effective rank is
measured on a FIXED held-out probe batch, which is the quantity that separated
D1's two arms (7.38 vs 1.76).

usage: decorr_causal.py <decorrelation_weight> <steps>
"""
import sys
import numpy as np, torch, torch.nn.functional as F
from torch import nn
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.preflight import restrict_cohort_to_split
from morpheus.v2.runner import attach_v2_targets, UncappedHoptimusBatches, _truncate_batch
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.training import V2LossSchedule, V2Trainer, PairedBiologyMemoryBank

# biology_full_consistency is hardcoded at 1.0 in the profile, and the G2.6
# isolation measured it blocking descent at that weight (1.847 vs 0.0034 with it
# off). It is therefore a second candidate driver and must be varied too.
_orig_weights = V2LossSchedule.weights


def _patched(self, epoch):
    w = _orig_weights(self, epoch)
    if w.get("biology_full_consistency"):
        w["biology_full_consistency"] = CONSIST
    return w


V2LossSchedule.weights = _patched

DECORR = float(sys.argv[1])
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 400
CONSIST = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
SEED = 42
D = "/home/ubuntu/e0_run/data/"
S = D + "paired_split_maximal.json"
torch.manual_seed(SEED); np.random.seed(SEED)
data = load_bio_query_data(D + "v1_abs_hallmark.json", S, wsi_mode="hoptimus_patch")
data, _ = restrict_cohort_to_split(data, S)
fit = np.asarray(data.split).astype(str) != "test"
attach_v2_targets(data, fit); data._v2_programme_head_dim = 256
split = np.asarray(data.split).astype(str)
train_idx, test_idx = np.where(split == "train")[0], np.where(split == "test")[0]

probe = _truncate_batch(next(iter(UncappedHoptimusBatches(data, test_idx, 8192, 7, shuffle=True))), 256)
probe = {k: (v.cuda() if isinstance(v, torch.Tensor) else v) for k, v in probe.items()}
NP = len(probe["indices"])

schedule = V2LossSchedule(objective_profile="programme_free", warmup_epochs=0,
                          variance_after_warmup=0.01, decorrelation_after_warmup=DECORR)
torch.manual_seed(SEED)
model = TumorStateV2(V2ModelConfig(hidden_dim=512, layers=4, heads=8), programme_dim=256).cuda()
opt = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-2)
trainer = V2Trainer(model, opt, schedule, "cuda", gradient_diagnostics_every=0)
trainer.biology_memory = PairedBiologyMemoryBank(capacity=4096)
trainer.prime_biology_memory(UncappedHoptimusBatches(data, train_idx, 8192, SEED, shuffle=False),
                             minimum_unique_keys=16)


def geometry():
    model.eval()
    with torch.no_grad():
        w = model(probe, view="wsi")["z_biology"].float()
        r = model(probe, view="rna")["z_biology"].float()
    model.train()
    wn, rn = F.normalize(w, dim=-1), F.normalize(r, dim=-1)
    eye = torch.eye(NP, dtype=torch.bool, device="cuda")
    # R3 (the historical probe statistic), reported beside canonical R1 so this
    # script's numbers can be read against both the older probes and the paper.
    return (effective_rank(w, variant=RANK_VARIANTS["R3"]), effective_rank(w),
            float(wn.std(0).mean()), float((rn @ rn.T)[~eye].mean()))


print(f"decorrelation={DECORR}  consistency={CONSIST}  steps={STEPS}  probe patients={NP}", flush=True)
print(f"{'step':>6}{'R3(hist)':>10}{'R1(canon)':>11}{'feat-std':>10}{'rna-rna':>9}{'contrastive':>13}{'decorr':>10}", flush=True)
step = 0
while step < STEPS:
    for raw in UncappedHoptimusBatches(data, train_idx, 8192, SEED + step, shuffle=True):
        batch = {k: (v.cuda(non_blocking=True) if isinstance(v, torch.Tensor) else v) for k, v in raw.items()}
        if step % 50 == 0:
            er, er1, sd, rr = geometry()
            print(f"{step:>6}{er:>10.2f}{er1:>11.2f}{sd:>10.4f}{rr:>9.4f}"
                  f"{last_c if step else float('nan'):>13.4f}"
                  f"{last_d if step else float('nan'):>10.2f}", flush=True)
        opt.zero_grad(set_to_none=True)
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            loss, metrics, _ = trainer.step(batch, epoch=1)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        last_c = metrics.get("biology_contrastive", float("nan"))
        last_d = metrics.get("decorrelation", float("nan"))
        step += 1
        if step >= STEPS:
            break
er, er1, sd, rr = geometry()
print(f"{step:>6}{er:>10.2f}{er1:>11.2f}{sd:>10.4f}{rr:>9.4f}{last_c:>13.4f}{last_d:>10.2f}", flush=True)
print(f"\nDONE decorrelation={DECORR} final_R3={er:.2f} final_R1_canonical={er1:.2f} "
      f"rna_rna={rr:.4f}", flush=True)
