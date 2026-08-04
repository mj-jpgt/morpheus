## 2026-08-04 10:00 UTC — The retraining envelope is 3.30× and bimodal; D1 becomes unresolvable; and the momentum separation is itself marginal against it

**Logged:** 2026-08-04 10:00 UTC. **How obtained:** five identical `programme_only` retrains
(`~/e0_run/d1_envelope/rep{1..5}`) read out by `v2/research/rebase/d1_envelope_readout.py`; six
seed-varied momentum runs (`~/e0_run/d1_diag/mseed_*`). Read against
`PREDECLARED_retraining_envelope_20260804T0330Z.md` (`bb36782`), committed before either reported.

### 1. The retraining envelope

Same seed 42, same configuration, five times. The only variation is GPU non-determinism.

| rep | rank (residualised) | channel |
|---|---:|---:|
| 1 | 28.320 | 0.6182 |
| **2** | **8.834** | **0.5859** |
| 3 | 28.348 | 0.6123 |
| 4 | 29.106 | 0.6110 |
| 5 | 28.959 | 0.6098 |

| | spread |
|---|---:|
| **rank** | **3.30×** |
| **channel** | **1.055×** |

**The distribution is bimodal, not a smooth spread.** Four runs agree to within 2% (28.32–29.11) and
one lands at a third of that (8.83). This is the same signature as the G2.6 gate's 6-of-8, and the
useful claim is not a range: **rank is reproducible ~80% of the time and catastrophically not ~20% of
the time**, on identical inputs. A range invites the reader to imagine a distribution; the shape says
what actually happens.

**Rep 2 is the cleanest single row this project has produced.** Same seed, same configuration: rank
falls 3.2×, channel falls 5%. One number moves by a factor; the other barely moves at all.

**The envelope is a floor, twice over**, as predeclared. `programme_only` is the stable arm — its
channel varies 1.018–1.026× across seeds where `programme_free` varies 1.056×, and its step-200
tripwire rank spans 1.18× across five seeds against `programme_free`'s 6.05×. And same-seed repeats
exclude seed variation entirely. A floor of 3.30× is a stronger statement than a point estimate of
3.30×.

### 2. Reading it against what was predeclared

The predeclaration fixed four bands. **3.30× falls in the `> 3.09×` band**, whose reading was written
before the number existed:

> No D1 rank difference is resolvable. D1 is **uninformative about rank**; the channel remains
> resolvable; the asymmetry is the finding. Explicitly **do not** claim the necessity result is
> refuted.

Applying it without amendment:

- **D1's rank ratios (2.02×, 3.09×, 1.68×) are all inside 3.30×.** D1 does not resolve whether rank
  tracks information, **in either direction**.
- **The necessity result is NOT refuted.** `programme_only` won the channel 3/3 with patient CIs
  excluding zero. That result stands at full prominence in §4.7, flagged as *not resolvable by this
  comparison* rather than deleted. A comparison inside the noise floor is not evidence for rank's
  reliability and equally not evidence against it.
- **The asymmetry is the finding.** On the *same five runs*, rank spreads 3.30× while the channel
  spreads 1.055×. On the *same three D1 seeds*, the channel differences have patient CIs excluding
  zero 3/3 while the rank differences do not clear the noise floor. The identical comparison shows the
  channel resolvable and rank not. That is the paper's thesis measured on one pair of arms — and it
  holds only because both sides are reported.

**§4.1's count moves from six of seven to seven of seven**, in our favour. Per the predeclaration this
gets the same scepticism as a result going the other way: it rests on **one** envelope measurement, of
**five** repeats, on **one** arm, at **one** seed, and it is a floor. It should be quoted as "seven of
seven fall inside a floor of 3.30× measured on five repeats", not as "seven of seven".

### 3. The momentum seed replication — and an awkward number

The seed-varied replication **did run** (`mseed_*`, 08:41, 500 steps, three seeds per momentum).
Values, canonical Roy & Vetterli order 1 and R3, held-out probe:

| m | seed 42 | seed 43 | seed 44 | spread |
|---|---:|---:|---:|---:|
| **0.999** canonical | 11.26 | 10.45 | 10.55 | 1.08× |
| **0** canonical | 3.18 | 1.13 | 2.36 | 2.81× |
| 0.999 (R3) | 7.40 | 6.85 | 7.15 | 1.08× |
| 0 (R3) | 2.81 | 1.05 | 2.06 | 2.68× |

**Every m=0.999 seed exceeds every m=0 seed on both statistics**, so §5.3's disjunction resolves in
favour of separation: the momentum fix is not an artefact of a single seed. The single-seed defect in
the original sweep is closed.

**But the separation is marginal against the envelope just measured.** Worst-case ratio, min(m=0.999)
over max(m=0): **10.45 / 3.18 = 3.29× canonical**, against a retraining floor of **3.30×**. On R3 it
is 6.85 / 2.81 = 2.44×, comfortably inside.

So by §4.1's own criterion, **the momentum fix's rank difference is not resolvable either.** This must
be stated rather than left for a referee.

Three caveats that make the comparison indicative rather than decisive, none of which rescue it:
the envelope was measured on `programme_only` at 40 epochs on the exported artifact; the momentum runs
are `programme_free` at 500 steps on a held-out probe. Different arm, duration and block. The
quantities are not substitutable, and a like-for-like envelope for this regime has not been measured.

**What does not depend on rank at all:** the momentum fix was adopted because the un-fixed
configuration collapses and the fixed one does not — a difference visible in retrieval, in the
contrastive loss, and in whether `programme_free` trains to 40 epochs without the tripwire firing. Its
justification does not rest on the rank ratio, and should not be written as if it does.

### In plain terms

Training the same model five times with everything identical gives the same rank four times and a
number three times smaller once. The information content barely moves at all across those same five
runs. That single comparison — one measure jumping by a factor while the other holds steady on
identical inputs — is the clearest statement of the paper's argument we have.

The cost is that our own experiment falls inside that noise. D1 cannot now tell us whether rank
tracks information; it can only tell us that the information difference is real and the rank
difference is not measurable. We do not get to say the inconvenient result went away — it is
unresolved, not refuted, and it stays in the paper at the size it was.

And the same standard bites our own fix: the momentum separation is 3.29× against a 3.30× floor. Not
the same arm or duration, so not a like-for-like disqualification — but close enough that pretending
otherwise would be exactly the double standard the paper accuses others of.

### Files / commits

- `~/e0_run/d1_envelope/rep{1..5}.npz`, `~/e0_run/d1_envelope_readout.log`
- `~/e0_run/d1_diag/mseed_m{0,0.999}_s{42,43,44}.log`
- `v2/research/rebase/d1_envelope_readout.py`, `d1_momentum_probe.py`
- Predeclaration: `PREDECLARED_retraining_envelope_20260804T0330Z.md` (`bb36782`)
