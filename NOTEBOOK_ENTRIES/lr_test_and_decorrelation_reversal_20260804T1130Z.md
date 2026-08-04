## 2026-08-04 11:30 UTC — The learning-rate test favours Account A but cannot rule out a third explanation my design omitted; and decorrelation's defect was conditional on the broken queue

**Logged:** 2026-08-04 11:30 UTC. **How obtained:** `~/e0_run/d1_diag/lr_L*.log` and
`ablate_decorr*.log`, 200 and 400 steps, `programme_free`, held-out probe. Read against
`PREDECLARED_learning_rate_test_20260804T2200Z.md` (`f68a7ac`), committed before the runs.

### 1. Learning-rate test — the predeclared reading, applied

| arm | lr | m | τ | eff-rank | rna_rna | predeclared A | predeclared B |
|---|---:|---:|---:|---:|---:|---|---|
| L1 | 1e-3 | 0.9 | 10 | **1.05** | 0.9257 | **fails** ✓ | works ✗ |
| L2 | 4e-5 | 0.99 | 100 | **27.88** | 0.5436 | **works** ✓ | fails ✗ |
| L3 | 1e-3 | 0 | 0 | 1.06 | 0.9946 | fails ✓ | fails ✓ |
| L4 | 4e-5 | 0.999 | 1000 | 35.24 | 0.3807 | works ✓ | works ✓ |

The predeclaration said: *"If L1 fails **and** L2 works, Account A."* L1 fails and L2 works, so on its
own terms **Account A (an absolute threshold in steps) is favoured and Account B (parameter-space
drift, τ_critical ∝ 1/lr) is falsified.** Both controls behaved as predicted, and **L4 — the
load-bearing control — passes**, so L2's success is not "the low-rate arm has not moved yet."

### But the design does not discriminate a third explanation, and I did not enumerate it

Learning rate separates the outcomes perfectly and momentum does not: both high-rate arms sit at
~1.05 regardless of m, both low-rate arms reach 27–35. Laying the design out as a 2×2 shows why:

| | m below threshold (0, 0.9) | m above threshold (0.99, 0.999) |
|---|---|---|
| **hi-lr 1e-3** | L1, L3 — **both fail** | *(not tested)* |
| **lo-lr 4e-5** | *(not tested)* | L2, L4 — **both work** |

**Every high-rate arm I ran has sub-threshold momentum, and every low-rate arm has supra-threshold
momentum.** The four results are therefore explained equally well by:

- **Account A** — m below threshold fails, above works, learning rate irrelevant; or
- **a pure learning-rate effect** — high rate collapses the representation regardless of m.

I designed L1/L2 as opposite-signed discriminators of A against B, and for that purpose they work. But
the diagonal cells are empty, so the design cannot separate A from an explanation I never wrote down.
**That is a defect in my design, not a property of the data**, and it would have been invisible
without the observation that lr separates the outcomes cleanly while m does not.

**The two missing cells are running** and are decisive:

| arm | lr | m | Account A predicts | pure-lr effect predicts |
|---|---:|---:|---|---|
| **L5** | 1e-3 | 0.999 | **works** (τ=1000 ≫ threshold) | **fails** |
| **L6** | 4e-5 | 0 | **fails** (τ=0 < threshold) | **works** |

Predeclared here, before they report. If L5 works and L6 fails, Account A survives a real test. If L5
fails and L6 works, the momentum threshold is an artefact of the learning rates it was measured at and
**the whole τ-threshold story goes the way of the other three**. Any mixed result falsifies both and
returns the mechanism to open.

Until L5/L6 report, the honest status is: **Account B is falsified; Account A is favoured but
confounded with learning rate; no mechanism is established.** This is the fourth account proposed for
this collapse and it has not yet earned more than the three that died.

### 2. Decorrelation — the earlier finding was conditional on the broken queue

All at m=0.999, 400 steps:

| decorrelation | eff-rank | rna_rna (mutual cosine) |
|---|---:|---:|
| 0.0 | 4.32 | 0.4774 |
| 0.01 | 6.22 | 0.7657 |
| 0.04 | **8.01** | **0.8696** |

**This reverses the earlier conclusion.** Without momentum, decorrelation *aggravated* collapse — 1.59
against 2.17 at m=0, with RNA cosine 0.98 against 0.42. With momentum it *raises* effective rank
monotonically, 4.32 → 8.01.

So *"`feature_decorrelation` is defective"* was never unconditional. It was **conditional on a queue
written by the query encoder**. Once the key encoder is decoupled, the same term helps by the rank
measure. The draft states the defect without that condition and must be corrected: every claim about
decorrelation needs "in the absence of a momentum key encoder" attached, and the reversal reported.

### And a dissociation stronger than anything currently in §4.9

Read the two columns together. As decorrelation rises, **effective rank rises 4.32 → 8.01** while
**the RNA-view mutual cosine also rises 0.477 → 0.870**. Rank says the representation is improving;
a direct collapse measure says it is degenerating — monotonically, in the same direction, on the same
three runs.

This is stronger than the existing §4.9 material for two reasons. It is **monotone across three
levels** rather than a single contrast, and the contradicting quantity is **co-measured on the
identical runs** rather than compared across tables. A mutual cosine of 0.870 means the RNA-view
states of different patients are nearly the same vector — the condition rank is supposed to detect —
and rank moves the wrong way as it worsens.

**Caveat, stated because the same standard applies here:** one seed per level, 400 steps. Given a
retraining floor of 3.30×, the 4.32 → 8.01 rank change (1.85×) is **inside the noise floor** and the
monotonicity across three levels is the load-bearing part, not the size. The cosine trend is the more
robust half of the dissociation.

### In plain terms

The learning-rate experiment killed the explanation it was designed to kill, and appeared to confirm
the survivor. But laying the runs out in a grid shows I only ever tested fast-learning with weak
momentum and slow-learning with strong momentum — so "momentum matters" and "learning rate matters"
fit the results equally well. The two missing combinations are running.

Separately, a penalty term we had written off as harmful turns out to be helpful once the queue is
fixed. And as we increase it, the rank score improves while a direct measurement of the thing rank is
supposed to detect gets steadily worse. That is the paper's argument in miniature, on three runs.

### Files / commits

- `~/e0_run/d1_diag/lr_L{1,2,3,4}_*.log`; L5/L6 in flight
- `~/e0_run/d1_diag/ablate_decorr{0.0,0.01,0.04}.log`
- Predeclaration: `PREDECLARED_learning_rate_test_20260804T2200Z.md` (`f68a7ac`)
- Earlier, now conditional: `d1b_premise_fails_all_five_arms_collapse`, `turnover_criterion_FALSIFIED`
