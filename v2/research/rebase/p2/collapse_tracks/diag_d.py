"""D1 G2.6 diagnostic D -- can the real encoder memorise 16 patients AT ALL?

diag_c: in-batch InfoNCE alone (no memory, no consistency, no regularisers)
plateaus at 2.3958 out of ln(16)=2.7726 and never moves again.  Before any
regularisation fix can matter, establish the ceiling of the clean task.

Sweeps learning rate and adds an explicit per-view variance floor to see
whether the plateau is a collapse basin or an optimisation failure.
Reports the quantities that decide it: positive cosine, worst negative
cosine, the InfoNCE margin actually required for <=0.10, and effective rank.
"""
import numpy as np, torch, torch.nn.functional as F
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.preflight import restrict_cohort_to_split
from morpheus.v2.runner import attach_v2_targets, UncappedHoptimusBatches, _truncate_batch
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.losses import variance_floor

D = "/home/ubuntu/e0_run/data/"; S = D + "paired_split_maximal.json"
SEED, N = 42, 16
torch.manual_seed(SEED)
data = load_bio_query_data(D + "v1_abs_hallmark.json", S, wsi_mode="hoptimus_patch")
data, _ = restrict_cohort_to_split(data, S)
fit = np.asarray(data.split).astype(str) != "test"
attach_v2_targets(data, fit); data._v2_programme_head_dim = 256
train = np.where(np.asarray(data.split).astype(str) == "train")[0]
it = iter(UncappedHoptimusBatches(data, train, 8192, SEED, shuffle=True))
prime, seen = [], set()
while len(seen) < 16:
    b = next(it); prime.append(b); seen.update(int(v) for v in b["indices"].tolist())
batch = _truncate_batch(next(it), N)
batch = {k: (v.cuda() if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
LAB = torch.arange(N, device="cuda")
EYE = torch.eye(N, dtype=torch.bool, device="cuda")
# margin against 15 in-batch negatives needed for InfoNCE <= 0.10
NEED = 0.07 * np.log(15 / (np.exp(0.10) - 1))
print(f"cosine margin required for in-batch InfoNCE <= 0.10 with 15 negatives: {NEED:.4f}")


def probe(model):
    model.eval()
    with torch.no_grad():
        w = model(batch, view="wsi")["z_biology"]; r = model(batch, view="rna")["z_biology"]
    model.train()
    wn, rn = F.normalize(w, dim=-1), F.normalize(r, dim=-1)
    cross = wn @ rn.T
    loss = 0.5 * (F.cross_entropy(cross / 0.07, LAB) + F.cross_entropy(cross.T / 0.07, LAB))
    pos = cross[EYE]
    worst = torch.stack([torch.cat([cross[i, :i], cross[i, i + 1:]]).max() for i in range(N)])
    sv = torch.linalg.svdvals(wn.float() - wn.float().mean(0))
    erank = float((sv.sum() ** 2) / (sv ** 2).sum())
    return dict(loss=float(loss), acc=float((cross.argmax(1) == LAB).float().mean()),
                pos=float(pos.mean()), worst=float(worst.mean()), margin=float((pos - worst).min()),
                ww=float((wn @ wn.T)[~EYE].mean()), erank=erank,
                std=float(wn.std(0).mean()))


def run(tag, lr, steps, var_weight=0.0, clip=None):
    torch.manual_seed(SEED)
    model = TumorStateV2(V2ModelConfig(hidden_dim=512, layers=4, heads=8), programme_dim=256).cuda()
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.0)
    print(f"\n--- {tag}  lr={lr} steps={steps} var_weight={var_weight} clip={clip}")
    marks = {0, 25, 50, 100, 200, 400, 800, 1500, 2500, steps}
    for step in range(steps + 1):
        if step in marks:
            p = probe(model)
            print(f"    {step:5d} loss {p['loss']:.4f} acc {p['acc']:.3f} pos {p['pos']:.4f} "
                  f"worst-neg {p['worst']:.4f} min-margin {p['margin']:+.4f} wsi-wsi {p['ww']:.4f} "
                  f"eff-rank {p['erank']:.2f} std {p['std']:.4f}")
        if step == steps:
            break
        opt.zero_grad(set_to_none=True)
        w = model(batch, view="wsi")["z_biology"]; r = model(batch, view="rna")["z_biology"]
        wn, rn = F.normalize(w, dim=-1), F.normalize(r, dim=-1)
        loss = 0.5 * (F.cross_entropy(wn @ rn.T / 0.07, LAB) + F.cross_entropy(rn @ wn.T / 0.07, LAB))
        if var_weight:
            t = w.shape[-1] ** -0.5
            loss = loss + var_weight * (variance_floor(w, target_std=t) + variance_floor(r, target_std=t))
        loss.backward()
        if clip:
            torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
        opt.step()


run("D1 clean in-batch InfoNCE, lr 1e-3", 1e-3, 2000)
run("D2 clean in-batch InfoNCE, lr 3e-3", 3e-3, 2000)
run("D3 clean in-batch InfoNCE, lr 1e-2", 1e-2, 2000)
run("D4 clean in-batch InfoNCE, lr 3e-3 + grad clip 1.0", 3e-3, 2000, clip=1.0)
run("D5 in-batch InfoNCE + per-view variance floor w=1.0, lr 3e-3", 3e-3, 2000, var_weight=1.0)
print("\nDIAG_D_DONE", flush=True)
