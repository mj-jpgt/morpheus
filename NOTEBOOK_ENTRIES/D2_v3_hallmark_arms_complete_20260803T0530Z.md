## 2026-08-03 05:30 UTC — D2 v3 Hallmark arms complete (3/3), PBS arms training

**Logged:** 2026-08-03 05:30 UTC. **How obtained:** Lambda A100 `150.136.45.194`,
`~/e0_run/d2_v3/logs/progress.log` (5-min sampler) and the `TRAIN_SUCCESS.json` files written by
`phase_d.run_d2` after its independent liveness assertion.

### Technical
All three Hallmark arms reached 40/40 epochs and passed the G2 liveness gate:

```
~/e0_run/d2_v3/d2_v3_s42/d2_h_seed42/TRAIN_SUCCESS.json
~/e0_run/d2_v3/d2_v3_s43/d2_h_seed43/TRAIN_SUCCESS.json
~/e0_run/d2_v3/d2_v3_s44/d2_h_seed44/TRAIN_SUCCESS.json
```

Timing, three concurrent streams: launched 01:13 UTC, epoch 0 at 01:28 (15 min preflight + the
800-step G2.6 overfit gate), 40/40 at **04:45 UTC — 3 h 32 min** for three arms in parallel.

Throughput tracked the co-tenant almost exactly, which confirms the earlier diagnosis was right:

| window | co-tenant loky workers | load avg | min/epoch |
|---|---:|---:|---:|
| 01:28–03:30 | 25–38 | 80–95 | ~5.6 |
| 03:30–04:20 | 9–18 | 43–50 | ~3.3 |
| 04:20–04:45 | 4 | 5–15 | **~2.2** |

A 2.5× swing driven entirely by another agent's joblib pools finishing. At the uncontended rate the
whole 6-arm sweep would have been ~3 h; contended it is tracking ~6 h.

PBS arms began their own liveness gates at 04:45, took ~25 min, and started epochs at ~05:10. At the
current ~3 min/epoch they should finish ~07:15 UTC, followed by six exports, CALIBRA G4 controls and
the paired bootstraps.

### In plain terms
Half the work is done and passed its health checks. The speed of the whole thing turned out to be a
straight function of how busy the shared machine was — when the other agent's jobs finished, ours
ran two and a half times faster with no change on our side.

### Meaning for the claim
Nothing yet — no arm has been scored. The Hallmark arms existing and passing liveness is the
precondition for the stratified readout, which is the deliverable.

### Files / commits
- `~/e0_run/d2_v3/d2_v3_s{42,43,44}/d2_h_seed*/` — checkpoints, metrics, liveness, TRAIN_SUCCESS
- `~/e0_run/d2_v3/logs/progress.log`
