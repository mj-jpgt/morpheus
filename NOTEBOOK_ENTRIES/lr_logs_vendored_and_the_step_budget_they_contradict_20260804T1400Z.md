## 2026-08-04 14:00 UTC — §5.2a's learning-rate logs are vendored; they agree with every number the section quotes and contradict its step budget

**Logged:** 2026-08-04 14:00 UTC. **How obtained:** `v2/research/rebase/p2/figures/extract_from_box.py`,
the single scripted path any figure datum may enter this repository by — one `tar` stream from
`ubuntu@150.136.45.194`, no per-file `scp`. **Not a re-run:** six files were copied and nothing was
recomputed.

### 1. What was closed

Draft §6.4 named this as the paper's second provenance gap and priced it at six file copies:

> `lr_L{1..6}.log` were not [vendored], so §5.2a's six rows are asserted against **the draft's own
> table** rather than against a log — the weakest of the three source kinds the audit supports, and
> the one this project has twice found insufficient. **This is a provenance gap in the paper's only
> established mechanism result.**

The six logs now sit in `v2/research/rebase/p2/figures/data/e0_run/d1_diag/` beside
`ablate_decorr*.log` and `mseed_*.log`, hashed into `data/MANIFEST.json` and covered by that
directory's `* -text` rule, which is what keeps the digests meaningful on a Windows checkout. Their
box filenames carry the arm's learning rate and momentum, so they are `lr_L1_hi_m0.9.log` …
`lr_L6_lo_m0.log`, not `lr_L{1..6}.log`.

| file | sha256 | bytes |
|---|---|---:|
| `lr_L1_hi_m0.9.log` | `034207a6154adcf61573e0fbbac4c5b6b10eb432cc7123380a1e1db4dbc825f9` | 711 |
| `lr_L2_lo_m0.99.log` | `56309100ab541d9cc09d692a365a62f3628bc35de872f337b9dbcf58708b6fc0` | 715 |
| `lr_L3_hi_m0.log` | `5c40ee6f29b98f91853c93cf6c5aa82396212630cb029ef70741d34e83305d7c` | 710 |
| `lr_L4_lo_m0.999.log` | `855bd75372631a5dd5dd643649c9273a3533d4a95b30533d0463f04f73d4e812` | 719 |
| `lr_L5_hi_m0.999.log` | `3e199aefd82f377a24a34519989356dc87a9c7919b64bfa276f6ba794e445b99` | 718 |
| `lr_L6_lo_m0.log` | `84019fc73cf92003d74806733be56a8ad712c5ef8ac95f841b7d8d40cffc6895` | 711 |

Verified byte-identical against `sha256sum` on the box. No other vendored file changed; the manifest
diff is these six entries and the `fetched_utc` stamp.

Audit rows 57–60 were then re-pointed from `kind: markdown` (a substring assertion against
`paper/P2_RANK_DRAFT.md` §5.2a) to `kind: probe_log`, read at the step the logs carry. **Every one of
the four values re-resolved to the recorded number exactly before the source was switched** — 12.30,
1.06, 35.24, 1.05 — so nothing moved; only the provenance did.

### 2. What the logs recover that no notebook entry had recorded

Each log's first line is the harness echoing its own resolved argv. That closes two holes the
previous entries left open: the six arms ran at **decorrelation 0.04** and **seed 42**. §5.2a is
therefore reproducible from this repository, which it was not before.

| arm | momentum | decorrelation | capacity | lr | seed | **steps** |
|---|---:|---:|---:|---:|---:|---:|
| L1 | 0.9 | 0.04 | 4096 | 1e-3 | 42 | **200** |
| L2 | 0.99 | 0.04 | 4096 | 4e-5 | 42 | **200** |
| L3 | 0 | 0.04 | 4096 | 1e-3 | 42 | **200** |
| L4 | 0.999 | 0.04 | 4096 | 4e-5 | 42 | **200** |
| L5 | 0.999 | 0.04 | 4096 | 1e-3 | 42 | **200** |
| L6 | 0 | 0.04 | 4096 | 4e-5 | 42 | **200** |

### 3. And what they contradict

**§5.2a says 400 steps. The logs say 200.**

The rank values agree exactly — L3 1.06, L1 1.05, L5 1.05, L6 12.30, L2 27.88, L4 35.24 are the
`R3-rank` column at step 200 of the six logs, and step 200 is the **last** row each log has, which is
what makes them the "final eff-rank" the section quotes. It is only the budget that disagrees, and it
disagrees in four places at once: §5.2a's provenance line, Appendix A's row, and both
`lr_test_and_decorrelation_reversal_20260804T1130Z.md` and
`learning_rate_is_the_mechanism_20260805T0100Z.md`.

Two things are worth saying about which side is likely right. The first is that the logs are the
primary artifact and the prose is a report of them. The second is that
`PREDECLARED_learning_rate_test_20260804T2200Z.md` fixed the reading rule in advance as *"centred
effective rank on the held-out probe at step 200"* — so the logs agree with the predeclaration and
the prose does not.

**Reported, not substituted.** This is recorded as a second entry in `floor_audit.json`'s
`known_source_disagreements`, which is the mechanism this project already has for exactly this: a
value that disagrees with its source is reported and never quietly replaced by either side. Nothing
in the audit quotes "400 steps". Correcting the prose belongs to §5.2a's owner and is not done here.

### 4. The same gap, one section earlier — §6.4 flags only `lr_L*`, and §5.2 was in the same position

Looking for the learning-rate logs turned up that **§5.2's own headline table was equally
un-vendored and not labelled so**. Twelve more files came down the same tar stream:

| files | what they are | audit rows they now carry |
|---|---|---|
| `long_m{0,0.9,0.99,0.999}.log` | §5.2's step-0-to-600 momentum table (launched at a **1,500**-step budget, read to 600) | 43, 46, 47 |
| `mom_{0,0.99,0.999}_d0.04.log` | §5.2 "measurement 2", the staleness falsification: rank **and** the key-to-encoder cosine on the same rows at step 100 | 48, 49 |
| `turn_cap{2048,4096,8192}_m0.9{,5,9}.log` | §5.2's turnover falsification, read at step 250 | 51, 52, 53 |

Eight rows moved from `kind: markdown` to `kind: probe_log`. **All sixteen values re-resolved to the
recorded number exactly** before any source was switched. It also settles the **reading step** for
five rows that recorded none, which is what a step-specific floor needs to be applied at all.

Two rows are deliberately left on `markdown`. Row 42 (§5.1's five-arm regulariser sweep) and row 50
(§5.2's capacity sweep, whose capacity-4,096 comparator reads 2.16 at step 150 in
`decorr_causal_0.04.log` — a log with a **six**-column header this audit's parser refuses by
design). Naming them is better than guessing at them.

### 5. And what the two families say about each other

**§5.2's table and §5.2's measurement 2 are different runs of the same configuration, and the
section presents them as one.** At step 100:

| m | `long_*` (the table) | `mom_*` (measurement 2) | **fold** |
|---:|---:|---:|---:|
| 0 | 1.62 | 2.58 | **1.593×** |
| 0.99 | 6.49 | 6.65 | 1.025× |
| 0.999 | 7.03 | 6.89 | 1.020× |

Same momentum, same decorrelation 0.04, same capacity 4,096, same default learning rate, the same
hardcoded seed 42 and the same step-0 state of 67.55. They differ in step budget and in GPU
non-determinism and in nothing else — which makes this an **n = 2 retraining spread on the fixed
held-out probe**, the block on which the paper had never measured one. On the m = 0 arm it is
**1.593×**.

Recorded in `known_source_disagreements`; neither family is substituted for the other, and §5.2's
prose is left to that section's owner. §5.2 also describes the sweep as "40 epochs = 583 steps"
where the logs say `steps=1500`.

### 6. Suite

426 passing, unchanged.
