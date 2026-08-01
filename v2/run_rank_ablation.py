"""Real-data biology-head rank ablation (paper T4).

Trains the MORPHEUS V2 biology head on the real held-out-cancer cohort in the
`programme_only` profile and measures the effective rank of ``z_biology`` on the
held-out test split, with the F-R2 feature-decorrelation term OFF (0.0) vs ON.
This is the on-data demonstration that F-R1 (normalize) + F-R2 (decorrelation)
lift the biology head out of collapse (baseline ~5-6 of 256).

Reuses the real data loader, batching, model and trainer — no re-implementation.
Run from the workspace where ``morpheus`` resolves to the fix clone and the data
symlinks resolve, e.g.:

    PYTHONPATH=$WS python -m morpheus.v2.run_rank_ablation \
        --data-config morpheus/configs/v1.json \
        --split-file /abs/tumor_state_heldout_cancer.json \
        --decorrelation-weight 0.04 --epochs 8 --token-budget 16384 \
        --out runs/rank_ablation/on.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch

from morpheus.src.training.train_bio_query_former import load_bio_query_data
from morpheus.v2.model import TumorStateV2, V2ModelConfig
from morpheus.v2.runner import UncappedHoptimusBatches, attach_v2_targets
from morpheus.v2.training import V2LossSchedule, V2Trainer


def effective_rank(x: torch.Tensor) -> float:
    x = x - x.mean(dim=0, keepdim=True)
    singular = torch.linalg.svdvals(x.double())
    singular = singular[singular > 1e-12]
    if singular.numel() == 0:
        return 0.0
    p = singular / singular.sum()
    return float(torch.exp(-(p * p.log()).sum()))


@torch.no_grad()
def collect_biology(indices: np.ndarray, data, model, device: str, token_budget: int, seed: int):
    """Return (wsi_biology reps [N,256], aligned global patient indices [N])."""
    model.eval()
    loader = UncappedHoptimusBatches(data, indices, token_budget, seed, shuffle=False)
    reps, idx = [], []
    for batch in loader:
        gid = batch["indices"].cpu().numpy()
        batch = {k: v.to(device) if isinstance(v, torch.Tensor) else v for k, v in batch.items()}
        # WSI biology view — the view the decorrelation term acts on and the one
        # the paper reports the collapse fingerprint for.
        reps.append(model(batch, "wsi")["z_biology"].float().cpu().numpy())
        idx.append(gid)
    model.train()
    return np.concatenate(reps, 0), np.concatenate(idx, 0)


def biology_rank_on(indices, data, model, device, token_budget, seed):
    reps, _ = collect_biology(indices, data, model, device, token_budget, seed)
    return effective_rank(torch.from_numpy(reps)), len(reps)


def molecular_prompting_probe(train_idx, test_idx, data, model, device, token_budget, seed, targets):
    """Ridge probe wsi_biology -> Hallmark targets; report pooled and within-cancer
    (macro-cancer) Pearson averaged over targets, on the held-out test split.

    This is the T4 secondary readout (ii): does recovering biology-head rank change
    the model's actual molecular-prompting specificity, or is rank decoupled from the
    (confounded) task score?

    `targets` is an (n_patients, K) matrix aligned to data.patient_ids, with NaN for
    patients lacking a target vector. The loader's own hallmark is train-fold-only and
    is a constant placeholder on the held-out-cancer test split, so callers must pass
    regenerated targets (frozen_rna_targets) that actually vary on test."""
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    from morpheus.v2.honest_metrics import macro_group_pearson, _pearson
    Xtr, itr = collect_biology(train_idx, data, model, device, token_budget, seed)
    Xte, ite = collect_biology(test_idx, data, model, device, token_budget, seed)
    cancers = np.asarray(data.cancers)
    Ytr_all, Yte_all = targets[itr], targets[ite]
    # Keep only patients (rows) with a full target vector; the target set covers a
    # subset of patients.
    mtr, mte = np.isfinite(Ytr_all).all(1), np.isfinite(Yte_all).all(1)
    Xtr, Ytr = Xtr[mtr], Ytr_all[mtr]
    Xte, Yte, cte = Xte[mte], Yte_all[mte], cancers[ite][mte]
    if len(Xtr) < 10 or len(Xte) < 10:
        return {"prompting_error": f"too few covered patients (train={len(Xtr)}, test={len(Xte)})"}
    # Standardize the (L2-normalized, small-magnitude) biology features before Ridge;
    # without this a large alpha shrinks predictions to a near-constant -> nan Pearson.
    scaler = StandardScaler().fit(Xtr)
    Xtr, Xte = scaler.transform(Xtr), scaler.transform(Xte)
    ridge = Ridge(alpha=10.0).fit(Xtr, Ytr)
    pred = ridge.predict(Xte)
    g = np.nanmean([_pearson(pred[:, j], Yte[:, j]) for j in range(Yte.shape[1])])
    w = np.nanmean([macro_group_pearson(pred[:, j], Yte[:, j], cte) for j in range(Yte.shape[1])])
    return {"prompting_pooled_pearson": float(g), "prompting_within_cancer_pearson": float(w),
            "n_targets": int(Yte.shape[1]), "n_train_covered": int(len(Xtr)), "n_test_covered": int(len(Xte))}


def load_aligned_targets(npz_path, data):
    """Build an (n_patients, K) target matrix aligned to data.patient_ids from a
    frozen targets npz (patient_ids + scores), NaN where a patient has no target."""
    t = np.load(npz_path, allow_pickle=True)
    pid_to_row = {str(p): i for i, p in enumerate(t["patient_ids"])}
    scores = np.asarray(t["scores"], dtype=np.float64)
    K = scores.shape[1]
    out = np.full((len(data.patient_ids), K), np.nan, dtype=np.float64)
    for i, pid in enumerate(data.patient_ids):
        r = pid_to_row.get(str(pid))
        if r is not None:
            out[i] = scores[r]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-config", default="morpheus/configs/v1.json")
    ap.add_argument("--split-file", required=True)
    ap.add_argument("--decorrelation-weight", type=float, default=0.04)
    ap.add_argument("--epochs", type=int, default=8)
    ap.add_argument("--token-budget", type=int, default=16384)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, default=0, help="subset patients per split for a fast smoke; 0 = full")
    ap.add_argument("--targets-npz", default="", help="frozen_rna_targets.npz (patient_ids+scores) for the T4 specificity probe; the loader's own hallmark lacks test coverage")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data = load_bio_query_data(args.data_config, args.split_file, wsi_mode="hoptimus_patch")
    attach_v2_targets(data)  # residual programme targets fit on the train fold only
    split = np.asarray(data.split)
    train = np.where(split == "train")[0]
    test = np.where(split == "test")[0]
    if args.limit:
        train, test = train[: args.limit], test[: args.limit]
    programme_dim = int(np.asarray(data.hallmark).shape[1])

    model = TumorStateV2(
        V2ModelConfig(hidden_dim=512, heads=8, layers=4, local_slots=32, slide_slots=16, patient_slots=8),
        programme_dim=programme_dim,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-2)
    schedule = V2LossSchedule(objective_profile="programme_only", warmup_epochs=1,
                              decorrelation_after_warmup=args.decorrelation_weight)
    trainer = V2Trainer(model, optimizer, schedule, device, gradient_diagnostics_every=0)

    record = {
        "decorrelation_weight": args.decorrelation_weight, "epochs": args.epochs,
        "token_budget": args.token_budget, "seed": args.seed, "programme_dim": programme_dim,
        "n_train": int(len(train)), "n_test": int(len(test)), "device": device,
        "biology_dim": 256, "history": [],
    }
    rank0, n = biology_rank_on(test, data, model, device, args.token_budget, args.seed)
    record["history"].append({"epoch": -1, "biology_effective_rank": rank0, "n_test": n})
    print(f"[init] biology effective_rank(test)={rank0:.2f} (n={n})", flush=True)

    for epoch in range(args.epochs):
        start = time.time()
        loader = UncappedHoptimusBatches(data, train, args.token_budget, args.seed + epoch)
        metrics = trainer.train_epoch(loader, epoch)
        rank, n = biology_rank_on(test, data, model, device, args.token_budget, args.seed)
        row = {"epoch": epoch, "biology_effective_rank": rank, "n_test": n,
               "train_programme": float(metrics.get("programme", float("nan"))),
               "train_decorrelation": float(metrics.get("decorrelation", float("nan"))),
               "seconds": round(time.time() - start, 1)}
        record["history"].append(row)
        print(f"[epoch {epoch}] biology effective_rank(test)={rank:.2f}  "
              f"programme={row['train_programme']:.4f}  decorr={row['train_decorrelation']:.4f}  "
              f"{row['seconds']}s", flush=True)
        out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(record, indent=2))

    record["final_biology_effective_rank"] = record["history"][-1]["biology_effective_rank"]
    # T4 secondary readout (ii): does rank recovery change molecular-prompting specificity?
    try:
        if not args.targets_npz:
            raise ValueError("pass --targets-npz frozen_rna_targets.npz (loader hallmark lacks test coverage)")
        targets = load_aligned_targets(args.targets_npz, data)
        probe = molecular_prompting_probe(train, test, data, model, device, args.token_budget, args.seed, targets)
        record["prompting"] = probe
        print(f"[probe] pooled={probe['prompting_pooled_pearson']:.4f}  "
              f"within_cancer={probe['prompting_within_cancer_pearson']:.4f}  "
              f"(n_targets={probe['n_targets']})", flush=True)
    except Exception as exc:  # keep the rank result even if the probe fails
        record["prompting_error"] = str(exc)
        print(f"[probe] FAILED: {exc}", flush=True)
    Path(args.out).write_text(json.dumps(record, indent=2))
    print(f"[done] final biology effective_rank={record['final_biology_effective_rank']:.2f} -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
