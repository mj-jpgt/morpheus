## Implementation Status & Synthetic Proof-of-Concept

The main method (F-R1 normalize + F-R2 covariance decorrelation) and the diagnostic
instrument (effective rank of the biology head) are **implemented and unit-tested** on the
`fix/biology-collapse` branch of the MORPHEUS repository. This section records what is
*done and verified locally* versus what remains *queued for the full Lambda run*, so no
claim in this paper is mistaken for a completed full-data result.

### Landed and CPU-verified (commits on `fix/biology-collapse`)

| Step | Commit | Change | Gate |
|---|---|---|---|
| **F-R0** | `6ac5e7e` | `v2/tests/test_stress_collapse.py`: `effective_rank()` (Roy–Vetterli) + 9 fast CPU tests (rank-collapse detection, KL/Gaussian-NLL finiteness, log-variance clamp non-saturation, no-op/grad-flow, anchor-gate saturation, slot accounting). | 9/9 pass; detector fails on a synthetic rank-2 matrix, passes on a full-rank one. |
| **F-R1** | `ab10be0` | `model.py`: normalize the biology state once so the Gaussian-NLL heads and the rank-spreaders (`z_biology`, KL, supcon) share one geometry (§ Method, T1). | stress suite green; `log_variance` stays off the +8 clamp. |
| **F-R2** | `13b6337` | `losses.feature_decorrelation` (Barlow-style off-diagonal **correlation** penalty) wired at weight $0.04$, gated in the `full` and `programme_only` profiles with a min-batch guard. | stress suite green (see below). |

### Synthetic proof-of-concept for T4 (unit-level, NOT the full-data claim)

On the synthetic collapse fixture in `test_stress_collapse.py` (a $B{=}24$ batch with
deliberately **low-rank programme targets**, `hidden=64`, `programme_only` profile, 60 steps),
the F-R2 decorrelation term raises the biology-head effective rank:

| Decorrelation weight $\lambda$ | Biology effective rank |
|---|---|
| $0.0$ (off) | $10.3$ |
| $0.04$ | $\mathbf{21.2}$ ($+10.9$) |
| $0.15$ | $21.3$ (saturated) |

This is a **synthetic, unit-scale** demonstration that the term is (a) not a no-op and (b) the
correct mechanism — it more than doubles rank where the low-rank target pull would otherwise
dominate, and saturates by $\lambda{=}0.15$, supporting $0.04$ as a conservative default. It is
**not** a claim about the seed-42 held-out cohort; the real biology rank is $\sim$5–6 of 256 and
the full-data recovery is the queued Lambda ablation (T4).

> **Engineering note (recorded for reproducibility):** the *first* F-R2 implementation used a
> plain covariance penalty and was a **silent no-op** — because `z_biology` is L2-normalized,
> its raw per-feature variance is $O(1/d)$ and the term contributed $\sim\!10^{-6}$ at
> $\lambda{=}0.04$ (rank delta $0.00$, parameter delta $4\times10^{-5}$). The stress test caught
> this; standardizing the features to unit variance (a scale-invariant *correlation* penalty,
> as in `whitened_cross_covariance`) fixed it. This is exactly the failure mode the effective-rank
> gate is designed to prevent, and a concrete instance of the paper's thesis that anti-collapse
> regularizers must be verified against a rank instrument, not assumed to work.

### Queued for the Lambda full run (not in this draft)

- **T4 multi-seed $\lambda$ ablation:** biology-head effective rank and control-adjusted
  within-cancer specificity, with vs. without F-R2, across seeds 42/43/44 on the real
  11-development / 21-held-out cohort — and whether rank recovery transfers to the benchmark
  or stays decoupled (as R4/C2 predicts).
- **F-R3 escalation** (RNA-paired biology InfoNCE with neighbour-KL → 0), only if F-R1+F-R2
  under-deliver on real data.
