## 2026-08-04 03:30 UTC — The turnover criterion is FALSIFIED: 4 of 5 predeclared predictions wrong, and the discriminating one inverted

**Logged:** 2026-08-04 03:30 UTC. **How obtained:** `~/ws_d1/momentum_test.py` on the A100
(`150.136.45.194`), five arms varying queue capacity and momentum, all from the same initialisation
(effective rank 67.55 at step 0), real streaming batches, lr 2e-4, decorrelation 0.04. Predictions and
reading rule were committed in `PREDECLARED_turnover_criterion_20260804T0130Z.md` (`01e38ea`)
**before the runs were launched**.

### Result

Reading rule, fixed in advance: works ≥ 5, fails ≤ 3.5, between = marginal.

| # | capacity | m | τ/T | s150 | s200 | s250 | **predicted** | **observed** |
|---|---:|---:|---:|---:|---:|---:|---|---|
| P1 | 2048 | 0.9 | 1.04 | 2.82 | 2.81 | 3.08 | marginal/works | **fails** ✗ |
| P2 | 2048 | 0.95 | 2.09 | 3.14 | 3.28 | 3.53 | works | **fails** ✗ |
| P3 | 4096 | 0.95 | 1.04 | 2.34 | 2.84 | 2.91 | marginal | **fails** ✗ |
| P4 | 8192 | 0.95 | **0.52** | 3.97 | 3.94 | 3.67 | **fails** | **best of the m=0.95 arms** ✗ |
| P5 | 8192 | 0.99 | 2.61 | 5.15 | 5.54 | 5.19 | works | works ✓ |

**Four of five wrong. P4, the discriminating prediction, is inverted** — the arm with the *lowest*
τ/T is the healthiest of the three m = 0.95 arms, where the criterion required it to be the worst.
Among those three arms, effective rank *increases* as τ/T *decreases*: 3.53 at τ/T = 2.09, 2.91 at
1.04, 3.67 at 0.52.

**Two of the three predeclared falsification conditions fired**: P4 did not fail, and outcomes track
`m` alone rather than the ratio. Per the predeclaration, the account is wrong.

### What actually organises the data

Momentum, essentially alone. Collecting every arm measured, at step 250–300:

| m | τ | capacity 2048 | 4096 | 8192 |
|---:|---:|---:|---:|---:|
| 0 | 0 | — | 2.29 | — |
| 0.9 | 10 | 3.08 | 1.85 | — |
| 0.95 | 20 | 3.53 | 2.91 | 3.67 |
| 0.99 | 100 | — | 5.77 | 5.19 |
| 0.999 | 1000 | — | 7.82 | — |

Rank is monotone in `m` and roughly flat in capacity once momentum is present. At fixed m = 0.95,
tripling the queue from 2048 to 8192 moves effective rank from 3.53 to 3.67 — no effect. At fixed
capacity 4096, raising m from 0.95 to 0.999 moves it from 2.91 to 7.82.

An **absolute** time-constant threshold sits somewhere between τ = 20 (fails) and τ = 100 (works),
and it sits in the same place at every capacity tested — 2048, 4096 and 8192. That is a real pattern
in the data and it is *not* a ratio.

**It is deliberately not being proposed as a mechanism.** This is the third explanation for this
collapse; regulariser weighting and MoCo staleness were both falsified by measurement, and the
turnover criterion has now been falsified by an experiment designed specifically to test it. Naming
"absolute τ ≈ 20–100 steps" as mechanism four on the strength of the same data that killed mechanism
three would be exactly the error the predeclaration exists to prevent. It is recorded as an
observation with a testable shape — vary the learning rate, which changes how far the encoder moves
per step without changing τ or T — and left there.

### What stands, and what does not

**Stands:** a momentum key encoder rescues the representation at fixed capacity, and the effect is
large (2.58 → 6.89 at step 100; 2.43 → 7.61 at step 500), monotone in m, and durable. **Does not
stand:** any account of *why*, including the turnover criterion I proposed and predeclared.

The earlier framing — that a query-written queue lets queries and keys co-move so collapse is free,
and a decoupled encoder removes that escape — survives this experiment untouched, because it makes no
prediction about *how much* decoupling is required. That is also its weakness: it is currently
unfalsifiable as stated, which is why it must not be written up as established.

### In plain terms

I predicted that the momentum needed would depend on how fast the queue empties, worked out five
specific predictions, wrote them down, and then ran the test. Four came out wrong, and the one
designed to be decisive came out backwards — the setting I said should fail worst did best of its
group.

What the data actually says is simpler and less satisfying: slower is better, and how big the queue is
barely matters once you have a slow encoder. There is a threshold, it looks like it is in the same
place regardless of queue size, and I do not know why. Writing down the wrong prediction in advance is
what makes it possible to say that cleanly instead of quietly re-fitting the story.

### Meaning for the claim

The fix is real, large and durable, and may be used. The mechanism is open — for the third time — and
`paper/QUEUE_ANCHORING.md` must say so rather than carrying a criterion that its own predeclared test
refuted. The dimensionless criterion is withdrawn.

### Files / commits

- Predeclaration: `01e38ea`, `NOTEBOOK_ENTRIES/PREDECLARED_turnover_criterion_20260804T0130Z.md`
- `~/e0_run/d1_diag/turn_cap{2048,4096,8192}_m*.log`, `long_m*.log`
- Requires revision: `paper/QUEUE_ANCHORING.md`
