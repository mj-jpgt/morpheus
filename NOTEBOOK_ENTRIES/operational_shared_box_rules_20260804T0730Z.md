## 2026-08-04 07:30 UTC — Operational rules: three errors of one class — a reading that looks like an answer and isn't

**Logged:** 2026-08-04 07:30 UTC. **How obtained:** two of my own mistakes on the A100
(`150.136.45.194`), all caught and all cheap, recorded because they are the same error three times and each
was avoidable given the one before.

### Technical

**1. Thread oversubscription (CPU).** Launching CPU analysis without thread caps on a shared 30-core
box lets every joblib worker spawn a full BLAS pool. These multiply against each other and against
concurrent training, and everything collapses into cache and scheduler contention. Measured: a test
that runs in **6.27 s** capped did not finish in **25 minutes** uncapped, ~250×. I misread that as an
expensive test and nearly proposed weakening it. The whole 275-test suite runs in **49 s** capped.

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
```

For GPU training jobs, 2 is a reasonable value rather than 1. `torch.get_num_threads()` defaults to
`nproc` = 30, so three concurrent runners ask for ~90 intra-op threads on 30 cores unless told
otherwise.

**2. GPU memory oversubscription.** Having learned rule 1, I checked CPU load before adding ten jobs —
load was 11/30, which looked fine — and did **not** check GPU memory. The card filled to 80/80 GB and
two runs died with `torch.OutOfMemoryError`. The jobs were ~5.6 GB each against ~1.8 GB for the
lighter ones, a difference I had not measured before launching.

**The rule, covering both:** *before adding jobs to a shared box, check every resource they contend
for, not the one that bit you last time.* Concretely, before launch:

```
cat /proc/loadavg                                          # CPU
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader   # GPU
```

and know the per-job footprint of what you are adding. A D1 training runner is ~6 GB; the momentum
probe is ~5.6 GB; a G2.6 gate check is ~1.8 GB.

**3. Priority under contention.** When the card is tight, the long-running experiment wins and the
diagnostic yields. On finding 78/80 GB used with ten diagnostics and three D1-B runners, I killed four
diagnostics rather than risk an OOM taking a training arm that had been explicitly preserved. A
diagnostic is minutes to re-run; a training arm is hours and, in this case, one of only two surviving
contrastive arms on the project.

**4. `pkill -f` matches your own command.** `pkill -f momentum_test.py` issued over ssh kills the ssh
session's own shell, because that shell's command line contains the pattern. This silently prevented a
relaunch once and killed a set of runs I meant to keep once. Match on the interpreter instead:

```
ps -eo pid,comm,args | awk '$2=="python" && /momentum_test/ {print $1}'
```

**5. A reading that looks like an answer and isn't — the most dangerous of the three.** Checking
whether the new in-run rank tripwire had fired, I queried `rank_tripwire_observed` and got an empty
list for every run. That reads as a confident negative: *the gate never triggered*. It was wrong. The
runner prefixes epoch metrics with `train_`, so the key is `train_rank_tripwire_observed`, and the
tripwire had fired and logged on every run that reached step 200.

Had I stopped there I would have reported an inert gate — the exact failure this project has spent two
days learning to detect. This is worse than the thread and memory errors above, because **an empty
result looks like data**. A missing file raises; a wrong key returns `[]`, and `[]` is
indistinguishable from "measured, found nothing" unless you check.

The rule: *when a query returns nothing, confirm the query can return something before believing the
nothing.* Concretely — list the keys that exist before asserting a key's absence means anything:

```python
print(sorted({k for row in rows for k in row}))   # what IS there
```

All three errors in this entry are the same class: a measurement that appeared to answer the question
and did not. Threads made a fast test look slow; unchecked GPU memory made a launch look safe; a
mistyped key made a working gate look dead.

### In plain terms

I made the same kind of mistake twice in a day: adding work to a shared machine after checking only
one of the things that machine can run out of. The first time it was processor threads, the second
time graphics memory. Both were cheap to fix and both were avoidable by asking "what does this job
consume, and how much of that is left?" rather than "does the machine look busy?".

### Meaning for the claim

Nothing scientific. It protects other agents' work on a shared box, and it protects this project's
own long-running jobs from being killed by its own diagnostics — which nearly happened, to an arm that
had been explicitly flagged for preservation.

### Files / commits

- `~/e0_run/d1_diag/probevar_*.log` (the OOM traces)
- Prior: `blas_oversubscription_not_slow_tests_20260803T1800Z.md`
