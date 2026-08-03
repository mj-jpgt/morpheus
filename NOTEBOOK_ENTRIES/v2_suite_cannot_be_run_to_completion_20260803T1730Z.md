## 2026-08-03 17:30 UTC — The v2 test suite cannot be run to completion: one test costs 5+ CPU-hours

**Logged:** 2026-08-03 17:30 UTC. **How obtained:** `pytest` on the A100 (`150.136.45.194`),
`/tmp/v2tests.log` and `/tmp/e0single.log`, plus `/proc/<pid>/stat` sampling.

### Technical

The brief states the suite is green at 168 v2 + 59 root tests. **The root suite reproduces (59
passed, 36 s). The full v2 suite does not complete.** Two attempts ran for over two hours each and
both stopped emitting after exactly 31 dots, i.e. inside test #32:

```
morpheus/v2/tests/test_e0_basis_transfer.py::test_zero_shared_biology_construction_does_not_pass
```

Run in isolation on an otherwise moderately loaded box it had still not finished after **25 minutes**.
It is **not deadlocked** — sampling `/proc/<pid>/stat` shows process state `R` and
`utime+stime` climbing by ~9,300 ticks per 5 wall-seconds, i.e. ~18.6 CPU-seconds per wall second
across BLAS threads. By that point it had consumed **~1.83M ticks ≈ 5 CPU-hours** in a single test
and was still going.

The cost is structural, not accidental: the test loops over four seeds, and each iteration builds a
module fixture and runs `_run_arms`, which is bootstrap- and permutation-heavy. Four full arm
computations per test, and the file contains several more tests of the same shape.

**It is not affected by anything changed in the G2.6 work.** `test_e0_basis_transfer.py` imports only
`morpheus.v2.calibra.e0_basis_transfer` (plus `calibra.gates`), and there is no `conftest.py` anywhere
in the tree, so there is no import path by which `losses.py`, `training.py`, `runner.py` or
`phase_d.py` could reach it.

**What was verified instead**, before launching D1:

| suite | result |
|---|---|
| `morpheus/tests` (root, 59) | **59 passed**, 36 s |
| `test_programme_free.py` + `test_phase_d_pairing.py` + `test_stress_collapse.py` | **31 passed**, 29 min |
| full `morpheus/v2/tests` (168) | **did not complete** — see above |

The three files above are the ones that import the changed modules, and include three new tests
pinning the centring behaviour, the small-batch guard, and the fact that G2.6 still grades the
*uncentred* number.

### In plain terms

The project's test suite is documented as passing, and the small fast half of it does. The larger
half contains a single statistical test that chews through five hours of processor time without
finishing, so in practice nobody can run the suite end to end and see it go green — which probably
explains why it is quoted as a number rather than as something routinely re-run.

Whether that test would eventually pass is unknown. It is definitely still running rather than
stuck, and it has nothing to do with the code changed here.

### Meaning for the claim

Nothing scientific. But the standing instruction "keep the full test suite green, run it before you
launch anything expensive" is **not currently satisfiable as written**, and any future agent that
takes it literally will burn hours discovering this. Either that test needs a reduced seed/bootstrap
count under a `slow` marker, or the instruction needs to name the subset that is actually runnable.
Recording it so the next person does not repeat the two-hour discovery twice, as this session did.

### Files / commits

- `v2/tests/test_e0_basis_transfer.py::test_zero_shared_biology_construction_does_not_pass`
- `/tmp/v2tests.log`, `/tmp/e0single.log` on the A100 (ephemeral; not persisted)
