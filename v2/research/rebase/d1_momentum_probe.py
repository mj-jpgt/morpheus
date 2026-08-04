"""Test (2): does a MOMENTUM key encoder rescue rank at FIXED capacity 4096?

The capacity sweep implicated the key set but confounds two things: how
INCONSISTENT the keys are, and how MANY negatives there are. This holds capacity
at 4096 -- same negative count -- and varies only whether keys are written by the
live query encoder or by an EMA copy of it. If rank recovers, inconsistency was
the cause. If it does not, negative count was.

It also MEASURES staleness directly instead of inferring it: at intervals, the
patients currently in the queue are re-encoded with the current query encoder and
the cosine between each stored key and its fresh re-encoding is reported as a
function of how many steps ago it was written. That is the decay curve. It makes
a falsifiable prediction: momentum should flatten it substantially at the same
capacity. If the curve barely moves and rank still recovers, the MoCo story is
wrong.

usage: momentum_test.py <momentum> <decorrelation> <steps>
       momentum = 0 means no momentum encoder (current behaviour)
"""
import sys, copy
import numpy as np, torch, torch.nn.functional as F
from torch import nn
from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.preflight import restrict_cohort_to_split
from morpheus.v2.runner import attach_v2_targets, UncappedHoptimusBatches, _truncate_batch
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.training import V2LossSchedule, V2Trainer, PairedBiologyMemoryBank

MOM = float(sys.argv[1])
DECORR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 300
CAPACITY = int(sys.argv[4]) if len(sys.argv) > 4 else 4096
LR = float(sys.argv[5]) if len(sys.argv) > 5 else 2e-4
SEED_ARG = int(sys.argv[6]) if len(sys.argv) > 6 else 42
SEED = SEED_ARG
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
opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-2)
trainer = V2Trainer(model, opt, schedule, "cuda", gradient_diagnostics_every=0)
trainer.biology_memory = PairedBiologyMemoryBank(capacity=CAPACITY)
trainer.prime_biology_memory(UncappedHoptimusBatches(data, train_idx, 8192, SEED, shuffle=False),
                             minimum_unique_keys=16)
# BOTH arms enqueue manually below, so the only difference between them is WHICH
# encoder writes the keys. Left to itself `step()` enqueues the query encoder's
# outputs internally, which would double-enqueue in the momentum arm and make the
# two arms structurally different rather than differing in one variable.
trainer.freeze_biology_memory = True

# Momentum key encoder: an EMA copy, deliberately NOT in the optimiser and with
# requires_grad off, so nothing can train it by accident.
key_encoder = None
if MOM > 0:
    key_encoder = copy.deepcopy(model).cuda()
    for p in key_encoder.parameters():
        p.requires_grad_(False)
    key_encoder.eval()

# The bank has no notion of when a slot was written, so track it alongside.
write_step = np.full(CAPACITY, -1, dtype=np.int64)


@torch.no_grad()
def ema_update():
    for pk, pq in zip(key_encoder.parameters(), model.parameters()):
        pk.mul_(MOM).add_(pq.detach(), alpha=1.0 - MOM)
    for bk, bq in zip(key_encoder.buffers(), model.buffers()):
        bk.copy_(bq)


@torch.no_grad()
def geometry():
    """Reports BOTH: the canonical Roy & Vetterli order-1 statistic that P2 quotes,
    and the R3 participation ratio the earlier sweep used, on the same states.
    Block: held-out test patients, WSI biology view, rows at own norms, centred."""
    model.eval()
    w = model(probe, view="wsi")["z_biology"].float()
    r = model(probe, view="rna")["z_biology"].float()
    model.train()
    wn, rn = F.normalize(w, dim=-1), F.normalize(r, dim=-1)
    eye = torch.eye(NP, dtype=torch.bool, device="cuda")
    # BOTH statistics come from the canonical function, neither is computed inline.
    # A fourth statistic was found living under the R2/R3 names in an analysis
    # script -- an inline order-2 Hill number of the EIGENVALUE distribution rather
    # than the repo's singular-value form. It changed a published count and survived
    # review, the tests and its own author. Importing is the only defence.
    from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank as _erank
    r3 = _erank(w, variant=RANK_VARIANTS["R3"])
    canonical = _erank(w)
    return (r3, float(wn.std(0).mean()), float((rn @ rn.T)[~eye].mean()), canonical)


@torch.no_grad()
def staleness(step: int, n_slots: int = 48):
    """cosine(stored key, fresh re-encoding of the same patient) vs slot age."""
    bank = trainer.biology_memory
    size = bank.size
    live = np.flatnonzero(write_step[:size] >= 0)
    if len(live) < 8:
        return []
    ages = step - write_step[live]
    order = live[np.argsort(-ages)]
    pick = order[np.linspace(0, len(order) - 1, min(n_slots, len(order))).astype(int)]
    rows = bank.indices[:size].cpu().numpy()[pick]
    keep = np.isin(rows, train_idx)
    pick, rows = pick[keep], rows[keep]
    if len(pick) < 8:
        return []
    try:
        batch = next(iter(UncappedHoptimusBatches(data, np.unique(rows), 8192, 3, shuffle=False)))
    except StopIteration:
        return []
    batch = {k: (v.cuda() if isinstance(v, torch.Tensor) else v) for k, v in batch.items()}
    model.eval()
    fresh = F.normalize(model(batch, view="wsi")["z_biology"].float(), dim=-1)
    model.train()
    got = batch["indices"].cpu().numpy()
    where = {int(v): i for i, v in enumerate(got)}
    out = []
    stored = F.normalize(bank.wsi_states[:size].float(), dim=-1)
    for slot, row in zip(pick, rows):
        if int(row) in where:
            out.append((int(step - write_step[slot]),
                        float(torch.dot(stored[slot], fresh[where[int(row)]]))))
    return out


tau = (1.0 / (1.0 - MOM)) if MOM > 0 else 0.0
print(f"momentum={MOM}  decorrelation={DECORR}  capacity={CAPACITY}  lr={LR}  seed={SEED}  steps={STEPS}  probe={NP}", flush=True)
print(f"tau={tau:.1f} steps   T=capacity/214={CAPACITY/214:.1f} steps   tau/T={tau/(CAPACITY/214):.2f}", flush=True)
print(f"{'step':>6}{'R3-rank':>10}{'CANONICAL':>11}{'feat-std':>10}{'rna-rna':>9}{'contrastive':>13}", flush=True)
step = 0
while step < STEPS:
    for raw in UncappedHoptimusBatches(data, train_idx, 8192, SEED + step, shuffle=True):
        batch = {k: (v.cuda(non_blocking=True) if isinstance(v, torch.Tensor) else v) for k, v in raw.items()}
        if step % 50 == 0:
            er, sd, rr, can = geometry()
            print(f"{step:>6}{er:>10.2f}{can:>11.2f}{sd:>10.4f}{rr:>9.4f}"
                  f"{(last_c if step else float('nan')):>13.4f}", flush=True)
            curve = staleness(step)
            if curve:
                buckets: dict[int, list[float]] = {}
                for age, cos in curve:
                    buckets.setdefault(min(age // 25 * 25, 250), []).append(cos)
                summary = "  ".join(f"{a}:{np.mean(v):.3f}" for a, v in sorted(buckets.items()))
                print(f"        STALENESS age:cos  {summary}", flush=True)
        opt.zero_grad(set_to_none=True)
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            loss, metrics, _ = trainer.step(batch, epoch=1)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        last_c = metrics.get("biology_contrastive", float("nan"))
        # Enqueue AFTER the optimiser step, from the key encoder. Both arms take
        # the same extra forward pass; only the ENCODER differs.
        if MOM > 0:
            ema_update()
            encoder = key_encoder
        else:
            encoder = model
            model.eval()
        with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            kw = encoder(batch, view="wsi")["z_biology"]
            kr = encoder(batch, view="rna")["z_biology"]
        if MOM <= 0:
            model.train()
        bank = trainer.biology_memory
        start = bank.cursor
        bank.update(kw.float(), kr.float(), batch["indices"])
        for j in range(len(batch["indices"])):
            write_step[(start + j) % CAPACITY] = step
        step += 1
        if step >= STEPS:
            break
er, sd, rr, can = geometry()
print(f"{step:>6}{er:>10.2f}{can:>11.2f}{sd:>10.4f}{rr:>9.4f}{last_c:>13.4f}", flush=True)
curve = staleness(step)
if curve:
    buckets = {}
    for age, cos in curve:
        buckets.setdefault(min(age // 25 * 25, 250), []).append(cos)
    print("        STALENESS age:cos  " + "  ".join(f"{a}:{np.mean(v):.3f}" for a, v in sorted(buckets.items())), flush=True)
print(f"\nDONE momentum={MOM} decorr={DECORR} final_eff_rank={er:.2f} rna_rna={rr:.4f}", flush=True)
