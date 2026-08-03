## 2026-08-03 18:30 UTC — The D1 pipeline's own output file is the contaminated readout, and it is the one that looks authoritative

**Logged:** 2026-08-03 18:30 UTC. **How obtained:** reading `run_d1` in
`v2/research/rebase/phase_d.py` against `d2_compare`'s `--target-groups` handling, while D1 trains
(epoch 15/40). Follows the 17:00 pre-registration entry; this is the mechanism by which that
pre-registration could be silently defeated.

### Technical

`run_d1`'s final stage invokes `d2_compare` like this:

```python
comparison = [sys.executable, "-m", "morpheus.v2.research.rebase.d2_compare",
              "--hallmark-artifacts", *exports["P"], "--pbs-artifacts", *exports["F"],
              "--targets", ..., "--output", str(root / "D1_PAIRED_BOOTSTRAP.json"),
              "--repeats", str(args.bootstrap_repeats),
              "--label-a", "programme_only", "--label-b", "programme_free", "--experiment", "D1"]
```

**No `--target-groups`.** `d2_compare` with that argument omitted scores *every non-control target* —
all 90 of `hallmark_in_training (50) + heldout_pathway (24) + immune_tme (8) + tumour_state (8)`.

So the pipeline will write `D1_PAIRED_BOOTSTRAP.json`, record it in `SUCCESS.json` under the key
`"paired_bootstrap"`, and that file will contain the readout in which **56% of the targets are
`programme_only`'s own training supervision**. It has the most authoritative-looking name in the run
root, it is the only bootstrap the pipeline produces, and it is the artifact a reader arrives at
first.

**This is precisely how D2's headline became uninterpretable** — the Stage-1 audit records that its
quoted PBS −0.10 result came from the equivalent unrestricted file, and that 50 of its 90 scored
targets were one arm's own supervision. The pre-registration I wrote at epoch 2/40 says the stratified
40-target readout is primary; that pre-registration is worth nothing if the pipeline emits only the
contaminated number under the obvious filename and someone later reads the run root without the
notebook.

**Consequences, and what is being done.**

1. The stratified and negative-control readouts are **not produced by the pipeline** and must be run
   afterwards by hand, from the same six exported artifacts:
   - `--target-groups heldout_pathway immune_tme tumour_state` → **the D1 result**
   - `--target-groups random_control` → the negative control both arms must fail
2. `D1_PAIRED_BOOTSTRAP.json` is to be read as the secondary, contaminated number, and quoted only
   with the 50-target caveat attached.
3. A `D1_READOUT_INDEX.json` will be written into the run root naming which file is the headline and
   why, so the run root is self-describing without the notebook.

**Not fixed in the running job, deliberately.** `phase_d.py` has been changed locally so future D1
runs emit all three readouts and record the stratified one as the headline, but that change is **not**
being synced to the box while the job is live. The running `phase_d` already holds its argv
construction in memory, so editing the file could not help this run anyway, and the only thing an
edit could do is introduce risk to a four-hour job for no benefit. It will be synced after D1
completes.

### In plain terms

The experiment's own machinery writes exactly one results file, gives it the most official-sounding
name in the folder, and lists it in the success record — and that file contains the version of the
result where half the exam questions were ones that one of the two candidates had been drilled on.

Deciding in advance which number counts, as I did earlier today, does not help if the machinery still
puts the wrong number in the obvious place. So the honest numbers get computed separately afterwards,
and a short index file goes in the folder saying which is which, so that someone opening the results
next month does not have to know that this note exists.

### Meaning for the claim

Nothing scientific yet. It removes the most likely route by which D1 would have repeated D2's
contamination — not through anyone deciding to quote the wrong number, but through the pipeline
making the wrong number the path of least resistance.

### Files / commits

- `v2/research/rebase/phase_d.py` — `run_d1`'s comparison stage (local change pending, not yet synced)
- `v2/research/rebase/d2_compare.py` — `--target-groups`, default "every non-control target"
- `~/e0_run/d1_v1/` — `D1_PAIRED_BOOTSTRAP.json` (will be written by the running job), `SUCCESS.json`
