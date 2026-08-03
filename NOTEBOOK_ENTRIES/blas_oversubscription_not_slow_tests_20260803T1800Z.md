## 2026-08-03 18:00 UTC — BLAS thread oversubscription, not an expensive test: a 6-second test took >25 minutes, and I misdiagnosed it

**Logged:** 2026-08-03 18:00 UTC. **How obtained:** `pytest` on the A100 (`150.136.45.194`), with and
without thread caps, plus `/proc/<pid>/stat` sampling and `/proc/<pid>/task` counts. **This entry
supersedes and corrects an earlier one of mine** (`v2_suite_cannot_be_run_to_completion`, commit
`0797e78`), which drew the wrong conclusion from a correct measurement.

### Technical

**The measurement.** Two full `v2/tests` runs stopped emitting after exactly 31 dots, inside
`test_e0_basis_transfer.py::test_zero_shared_biology_construction_does_not_pass`. Run in isolation it
had still not finished after 25 minutes. Sampling `/proc/<pid>/stat` showed process state `R` and
`utime+stime` climbing ~18.6 CPU-seconds per wall second — ~1.83M ticks (~5 CPU-hours) consumed and
still going.

**My conclusion was wrong.** I read "5 CPU-hours burned" as "this test is inherently expensive",
concluded the suite could not be run to completion, and proposed a `slow` marker or a reduced
seed/bootstrap count. All of that was wrong, and proposing to weaken a test on the strength of a
misread environment is the more serious of the two errors.

**The actual cause.** ~19 threads running flat out is not a heavy test; it is oversubscription
thrashing. The correct reading of the same number: the box was at load 110–130 across **30 cores**,
carrying concurrent D1 training and another agent's CALIBRA work. joblib workers each spawn a full
BLAS thread pool, those multiply against every other pool on the box, and the whole thing collapses
into cache and scheduler contention. Bootstrap SVDs — exactly what this test does — degrade worst
under that.

**Verified on the same contended box, same test:**

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  ~/venv/bin/python -m pytest ...::test_zero_shared_biology_construction_does_not_pass -q
1 passed in 6.27s
```

**6.27 seconds** capped, against not finishing in 25+ minutes uncapped. Roughly a 250× swing from
thread capping alone, on the same machine at the same moment. Independently, the coordinator measured
the same test at 10.33 s and the whole 275-test suite at 58 s on an uncontended machine.

**The suite is fine and remains green — now 275 tests**, up from 168 with the new tests from this
session and from the other agents. No `slow` marker, no reduced seed count, and no change to the
standing instruction is warranted; the earlier entry proposing them is withdrawn.

**The same pathology is on the training side.** Each D1 runner holds **123 threads** with no caps in
its environment, and `torch.get_num_threads()` defaults to 30 = `nproc`. Three concurrent runners
therefore ask for ~90 intra-op threads on 30 cores, and the export / CALIBRA / paired-bootstrap
stages that `run_d1` launches later inherit the same uncapped environment — CALIBRA at
`--calibra-jobs 6`, each worker spawning its own 30-thread pool, is the textbook case. **This is
plausibly part of why D2's paired bootstrap was measured at ~96 CPU-hours.** D1 was *not* restarted to
fix it: it is 8 epochs in, the runs are correct either way, restarting costs the progress and the
`phase_d` bookkeeping the audit depends on, and box load has since fallen to ~19 as the other agent's
work finished. Worth capping deliberately on the next launch rather than retrofitting this one.

**Standing practice from here:** every CPU analysis launched on a shared box gets
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`, and any future
`phase_d` invocation should carry them in its environment so the analysis stages inherit them.

### In plain terms

A test that takes six seconds appeared to take over twenty-five minutes, and I concluded the test was
enormous and should be trimmed. It was not. Every program on the machine had been told it could use
all thirty processors at once, so a dozen of them were fighting over the same thirty processors and
all of them crawled. Telling the test to use one processor made it finish in six seconds.

The lesson worth keeping is the one about method, not about threads: I had a correct measurement and
reached for the explanation that blamed the code, when the explanation was the environment I was
measuring in. The tell was there — nineteen threads at full tilt is a symptom of contention, not of a
big computation — and I read it the wrong way round.

### Meaning for the claim

Nothing scientific. Two corrections to the record: the suite is green at 275 and always was, and my
proposal to weaken a regression test — one whose docstring says it exists because the pre-fix
implementation could not produce the number zero — is withdrawn before it was acted on. Nothing was
changed in that test.

### Files / commits

- Supersedes `0797e78` (entry deleted in this commit; the erroneous text remains in git history)
- `v2/tests/test_e0_basis_transfer.py` — **unchanged**, no marker added, no counts reduced
