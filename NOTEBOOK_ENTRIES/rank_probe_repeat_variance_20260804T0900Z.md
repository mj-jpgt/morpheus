## 2026-08-04 09:00 UTC — The rank probe is tight and cleanly separated, but 3 repeats cannot rule out the gate's failure mode

**Logged:** 2026-08-04 09:00 UTC. **How obtained:** `~/ws_d1/momentum_test.py` on the A100
(`150.136.45.194`), 200 steps, real streaming batches, live queue, **identical inputs**, repeated per
condition. Logs `~/e0_run/d1_diag/probevar_*.log`.

### Technical

Centred effective rank on the held-out probe, three repetitions per condition:

| step | m = 0.999 (rep 1/2/3) | m = 0 (rep 1/2/3) |
|---:|---|---|
| 0 | 67.55 / 67.55 / 67.55 | 67.55 / 67.55 / 67.55 |
| 50 | 9.34 / 9.31 / 9.33 | 4.10 / 4.11 / 4.13 |
| 100 | 6.98 / 6.62 / 7.01 | 2.22 / 3.20 / 3.12 |
| 150 | 7.10 / 6.61 / 6.32 | 1.74 / 1.80 / 2.40 |
| **200** | **7.15 / 6.92 / 7.25** | **1.80 / 1.46 / 1.98** |

| | m = 0.999 | m = 0 |
|---|---|---|
| range at step 200 | 6.92 – 7.25 | 1.46 – 1.98 |
| relative spread | **4.7%** | 30% |
| **separation** | **3.5×, with an empty band from 1.98 to 6.92** | |

**Compared with G2.6 on the same stack:**

| | G2.6 | rank probe |
|---|---|---|
| repeats | 8 | 3 per condition |
| range | **650×** | 1.05× (working arm) |
| shape | **bimodal** — 6 tight, 2 divergent | unimodal, no outliers |
| worst outlier | 5.585, twice chance | none |
| steps | 2,400 | 200 |

The probe is dramatically better behaved on this evidence: it separates the two conditions by 3.5×
with a wide empty band, and a threshold at 4 sits in the middle of that band rather than near either
distribution.

### The caveat, and it is not small

**Three repetitions cannot rule out the gate's failure mode.** G2.6's divergence rate is 2/8 = 25%.
Observing zero divergences in three probe runs is entirely consistent with that rate — P(0 in 3 | p =
0.25) = 0.42 — and the exact upper 95% bound from 0/3 is **p ≤ 0.63**. So this experiment constrains
the *typical* spread well and the *tail* essentially not at all.

The planned design was ten repetitions. It was cut to six because the ten-way launch exhausted GPU
memory and killed two runs; four were then stopped deliberately to protect D1-B's surviving arms from
an OOM. That was the right operational call and it is a real weakening of this result, recorded rather
than glossed.

**What follows regardless of the tail.** The probe is 12× cheaper than the gate (200 steps versus
2,400), so repeat-and-minimise costs little on it. Whatever the true divergence rate, reading the
probe as best-of-N is cheap insurance and is the same discipline being adopted for G2.6. A replacement
gate should not be trusted on a tighter evidential basis than the one it replaces — and on divergence
rate specifically, it currently *is* on a thinner basis, not a better one.

### In plain terms

The proposed replacement check gives nearly the same answer every time it is run — the healthy setting
lands between 6.9 and 7.3, the broken one between 1.5 and 2.0, and nothing came anywhere near the line
between them. The old check, by contrast, ranged over a factor of six hundred and occasionally
declared a healthy model dead.

The honest limit: I ran the new check three times per setting, not ten, because the ten-way launch ran
the graphics card out of memory and I stopped several rather than risk killing a training run that had
been set aside for preservation. Three repeats show the usual spread is small; they say almost nothing
about how often it might go badly wrong, which is exactly the failure the old check has. So the new
check should be read as best-of-N too, not trusted because it looked tidy three times.

### Meaning for the claim

Both preconditions requested before implementation are now measured: G2.6's repeat variance
(650×, bimodal, 25% divergence) and the probe's (tight, unimodal, cleanly separated, tail
unconstrained). The evidence supports both proposed changes — repeat-and-minimise as the general
estimator rule, and the training-scale probe for D1 — with the probe's tail explicitly flagged as
under-measured.

### Files / commits

- `~/e0_run/d1_diag/probevar_m{0.999,0}_{1,2,3}.log`
- `~/e0_run/d1_diag/gatevar_1..8.log`
- Prior: `g26_is_not_reproducible_20260804T0700Z.md`, `operational_shared_box_rules_20260804T0730Z.md`
