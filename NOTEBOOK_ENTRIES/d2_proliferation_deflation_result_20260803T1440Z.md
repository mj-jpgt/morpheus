# Proliferation-deflated D2 — the gap does not shrink, it grows to 117–199% of baseline. But the rank-matched placebo reproduces that almost exactly, and the negative control degrades 2–6×, so this run is a clean NULL FOR PROLIFERATION and not a clean "gap survives" certificate. D3 now runs on all three seeds: 92.6–98.9% retained.

**Logged:** 2026-08-03 14:40 UTC. Pre-registered in
`NOTEBOOK_ENTRIES/d2_proliferation_deflation_preregistration_20260803T1400Z.md`, committed `3e6333f`,
**before** the deflated comparison ran. Code in `3d3f272`.

**How obtained:** Lambda box `150.136.45.194`, `~/ws_d3`, CPU only, thread-capped, GPU untouched (D1
training at 99% throughout). 30 `d2_compare` invocations, 2,000 paired patient+cancer bootstrap
repeats each, one seed-pair per invocation with `--seed 42` so the bootstrap RNG matches the existing
ledger runs exactly. Outputs `~/e0_run/d3/d2_deflate/bootstrap/`, axis cuts in `.../axes/`.

---

### Technical

#### Provenance check first: the baseline reproduces the ledger exactly

Re-running the undeflated comparison through the new code reproduces `D2_RESULT.md` to four decimals:
**−0.1325 / −0.1089 / −0.1226** on the untrained-40 targets and **−0.0099 / −0.0280 / −0.0268** on
the random controls, with the same CIs. So any movement below is the deflation, not a code change.

#### The cut, and which side it was applied to

The PBS **axes** were cut, not the readout target block. The proliferation axes' patient scores were
appended to the confound design, so `cross_fitted_residuals` removed that subspace from the targets
and from **both** arms identically before any canonical correlation. The readout stays the
pre-registered 40 `heldout_pathway + immune_tme + tumour_state` targets. `PBS_001` — the
4.5×-background axis that was the most legible axis in 4/6 D2.3 runs — is in all three real cuts.

Design goes 108 → 141 columns (33 axes), 140 (32), 154 (46).

#### Primary readout — `pbs − hallmark`, untrained-40 targets

| cut | n axes | seed 42 | seed 43 | seed 44 | patient CI₉₅ excludes 0 |
|---|---:|---:|---:|---:|---|
| baseline | 0 | −0.1325 | −0.1089 | −0.1226 | 3/3 |
| **`prol_top100`** (governing) | 33 | **−0.1962** | **−0.2164** | **−0.1844** | 3/3 |
| `prol_wmean` (co-primary) | 32 | −0.1850 | −0.2017 | −0.1440 | 3/3 |
| `union` (conservative) | 46 | −0.1966 | −0.2019 | −0.1482 | 3/3 |
| **`placebo_random`** (rank-matched) | 33 | **−0.1931** | **−0.1928** | **−0.1848** | 3/3 |

Retention of |Δ| against the same seed's baseline:

| seed | `prol_top100` | `prol_wmean` | `union` | **placebo** |
|---|---:|---:|---:|---:|
| 42 | 148% | 140% | 148% | **146%** |
| 43 | 199% | 185% | 185% | **177%** |
| 44 | 150% | 117% | 121% | **151%** |

**On the pre-declared primary bar the gap survives outright** — ≥70% retained in 3/3 seeds under
every cut including the conservative `union`, with the patient CI₉₅ excluding zero in 3/3 everywhere.
The gap does not shrink at all; it roughly doubles.

#### …and the placebo says the doubling is not proliferation

**The rank-matched placebo — 33 axes drawn at random from the 82 in neither proliferation quartile —
reproduces the widening almost exactly: 146% / 177% / 151% against the proliferation cut's 148% /
199% / 150%.** On seed 44 the placebo widens the gap *more* than the real cut does.

The per-arm decomposition says why. Removing any 33 PBS axes costs the PBS arm 1.6–2.3× more than the
Hallmark arm — unsurprising, because the PBS arm was *supervised* on those axes, so its channel lives
disproportionately in that subspace. Proliferation has nothing to do with it: the placebo ratios
(1.74–2.27×) sit inside the proliferation ratios (1.56–2.29×).

**So the widening must not be reported as "removing proliferation makes PBS look even worse."** That
is precisely the claim the placebo exists to prevent, and it is false.

#### The negative control degrades, and that is the reason this is not a certificate

Pre-declared: *"If deflation inflates the random-control gap towards the real one, the instrument has
been damaged by the deflation and no verdict is issued."* It does.

| cut | random-control Δ (42/43/44) | ratio real ÷ control |
|---|---|---|
| baseline | −0.0099 / −0.0280 / −0.0268 | **13.4× / 3.9× / 4.6×** |
| `prol_top100` | −0.0647 / −0.0878 / −0.0625 | 3.0× / 2.5× / 3.0× |
| `placebo_random` | −0.0654 / −0.0533 / −0.0503 | 3.0× / 3.6× / 3.7× |

The arm gap on 90 **random** gene sets inflates 2–6×, and its patient CI₉₅ — which included zero in
2/3 seeds at baseline — now excludes zero in 3/3. The instrument's discrimination between real and
random targets falls from 13.4/3.9/4.6× to about 3×. The placebo shows this damage is *also* generic.

A deflated comparison whose negative control has itself moved that far is **not the same instrument**
as the baseline one, and a "the gap survived" certificate read off it would be over-claiming.

#### What this run does establish, cleanly

**The proliferation cut is statistically indistinguishable from a rank-matched arbitrary cut on every
quantity measured** — gap retention, per-arm cost, negative-control inflation. That is a
placebo-controlled null: **there is nothing special about the proliferation programme in the D2 gap.**
It is a stronger and more honest statement than "the gap survives", because it does not depend on the
deflated instrument being undamaged.

#### D3 — now on all three seeds

Same command as seed 42, ABSOLUTE purity, complete-case n = 2,650, before/after on identical patients,
2,000 within-cancer permutations. `gates_pass: true` and `rna_positive_control_passed: true` in 3/3.

`wsi_biology`, `excess_over_null_median` retained (pre-declared bar ≥80%):

| seed | d2_h | d2_i |
|---|---:|---:|
| 42 | 94.2% | 97.9% |
| 43 | 96.3% | 98.9% |
| 44 | 92.6% | 98.0% |

**6/6 cells clear the bar**, observed still above `null_p95` in 6/6, `permutation_p` at the 1/2001
floor in 6/6. **The falsifier did not fire in any seed.**

The comparison that actually kills the tumour-content hypothesis, now on three seeds — morphology vs
the RNA→RNA circular control, same run, same patients:

| seed | d2_h: wsi loses / rna→rna loses | d2_i: wsi loses / rna→rna loses |
|---|---|---|
| 42 | 5.8% / 6.8% ✓ | 2.1% / 3.9% ✓ |
| 43 | 3.7% / 5.0% ✓ | 1.1% / 4.1% ✓ |
| 44 | **7.4% / 6.4% ✗** | 2.0% / 4.6% ✓ |

**The image channel loses less than the circular control in 5 of 6 cells, not 6 of 6.** The exception
is d2_h seed 44. The direction is consistent and the d2_i arm is unanimous, but the claim must be
quoted as 5/6 — purity shaves a common few percent off everything rather than selectively deflating
the image result, and one cell runs the other way.

### In plain terms

The worry was that our headline D2 result — that supervising on the perturbation basis works *worse*
than supervising on curated pathways — might just be a story about cell division. So we removed the
proliferation part of the perturbation basis from both arms and from the targets, and looked again.

The gap did not shrink. It roughly doubled. That looks like a resounding answer until you run the
control we built for exactly this: remove an equal-sized *random* chunk of the basis instead. That
doubles the gap too, by the same amount. The reason is mundane — one of the two arms was trained on
those axes, so taking any of them away hurts that arm more, whichever ones you take.

So the honest result is a negative about proliferation rather than a positive about anything:
proliferation is not doing anything here that an arbitrary equal-sized chunk of the basis doesn't also
do. And we cannot hand out a "the gap survived" certificate from this run, because the deflation also
degraded our sanity check — the difference between the arms on *random* gene sets grew several-fold,
which means the deflated measurement is a blunter instrument than the one it is being compared to.

### Meaning for the claim

**I am still not flipping `proliferation_controlled` for E0's `transfer` claim**, and
`tests/test_claim_guards.py` is untouched. Two reasons, the first of which was written down *before*
these numbers existed:

1. **Pre-declared direction argument** (pre-registration §5). D2's finding is that PBS supervision is
   *worse* than Hallmark. Showing a *negative* result about the perturbation basis is not
   proliferation-driven does not establish that E0's *positive* alignment is more than proliferation.
   Opposite signs, different claims. I committed this reading in advance precisely so it could not be
   fitted to whichever answer arrived.
2. **New, data-driven reason.** The deflated instrument's negative control degrades 2–6×, so this run
   cannot certify anything on its own — a point my own pre-registration named as a no-verdict
   condition.

**What IS discharged, and should be recorded as such:** the proliferation confound for **the D2 gap**
and for **PBS-axis legibility**, on two independent placebo-controlled results —
- this run: the proliferation cut is indistinguishable from a rank-matched random cut, and
- D2.3 (`d2_3_per_axis_proliferation_20260803T1345Z.md`): 85–95 of the 95 non-proliferation axes stay
  legible from morphology at ~90% of median.

**The one remaining step that would genuinely discharge the blocker for E0**, stated concretely so
nobody has to rediscover it: re-run E0's own statistic with **the responsive arm stratified by
proliferation loading** — `claim_guards`' remedy #2 — and show the responsive-minus-control gap
(+0.0727 at k=10, `E0_RESULT.md` §0) survives when proliferation genes are removed from the
responsive arm. That is a re-run of `runs/e0_20260731` (901 s on an A100, 1,000 Haar draws × 1,000
bootstraps, k ∈ {10,25,50,100}), it needs the GPU, and it was not in this task's scope. Until it is
run, `proliferation_deflation` stands for E0.

**Structural note worth acting on.** Nothing in production builds an E0 claim dict. The project's
record of E0's admissibility exists **only** as a hardcoded fixture at `tests/test_claim_guards.py:135`,
and `validate_claim` reads real evidence from nowhere. So "discharging a blocker" currently means
editing a test fixture by hand — which is fragile in both directions, and is worth replacing with an
evidence file that the guard actually reads.

### Limitations

1. The deflation is not arm-neutral: it removes a subspace one arm was supervised on. The placebo
   quantifies this but does not fix it. A genuinely neutral test would deflate a subspace neither arm
   was trained on.
2. `prol_wmean` and `union` on seed 44 retain 117%/121%, materially lower than the other seeds and
   than that seed's own placebo (151%). Seed 44 is the least stable cell; it still clears every bar.
3. Only the patient bootstrap is used for the pass/fail call. The cancer-cluster CI excludes zero in
   3/3 for every real cut on untrained40, and is reported, but it is the wider interval.
4. D3 remains seed-42-code for the placebo arm only (seeds 43/44 have no placebo run); the seed-42
   placebo is the rank control for all three.

### Files / commits

- `~/e0_run/d3/d2_deflate/{bootstrap/*.json,axes/*.txt,axes/axis_sets.json,jobs.txt,one.sh,logs/}`
- `~/e0_run/d3/main_absolute_seed4{3,4}/` (D3 seeds 43 and 44)
- Code `3d3f272`; pre-registration `3e6333f`.
