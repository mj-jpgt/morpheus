## 2026-08-03 23:50 UTC — My GPU workspace had silently drifted from HEAD, and it broke the negative control

**Logged:** 2026-08-03 23:50 UTC. **How obtained:** the D1-B audit failed on `random_control` with
*"target-group selection ['random_control'] left no targets"* — a bug whose fix was already committed
and present in my local checkout. Diagnosed by comparing checksums between the repository and the box.

### Technical

The GPU workspace `~/morpheus-rebase-d1` (reached as `~/ws_d1/morpheus`) was created once from a
tarball and has since been maintained **file by file**: whenever I changed something I `scp`'d that
file. Other agents' commits never reached it. Over two days it accumulated an unknown number of stale
files, and I did not notice because everything I personally touched was current.

It surfaced as a wrong scientific result, not as an error. `d2_compare._targets` had a documented fix
— when `--target-groups` is given, select by group only, because applying the `RANDOM_CONTROL__`
prefix drop on top of an explicit group request made `--target-groups random_control` select nothing.
That fix is in the repository. It was **not** in my workspace:

```
local  062470cb5daec1f449799079d2bb97ff  v2/research/rebase/d2_compare.py
box    03021de3f9005d391721d2773ca5785a
grep "Selection is BY GROUP ONLY" on the box -> 0 matches
```

So the negative control — **the one readout that proves the instrument is not manufacturing a channel
out of noise** — could not run, and the audit exited non-zero for a reason that had nothing to do with
the data.

**Why the earlier partial syncs did not catch it.** After the rank consolidation I synced the 15 files
in `git diff --name-only <my-last-commit>..HEAD -- v2/`. That is the set *changed since my last commit*,
which is not the set *differing from my workspace*. `d2_compare.py` had been fixed **before** that
range and so appeared in neither, while being stale on the box the whole time.

**Fixed** by shipping every tracked file (296) rather than a diff, and committing the workspace's
local-only git repo so drift is visible next time as a dirty tree. Verified: `clean=0`, and the fix
now greps.

**One cosmetic side effect, recorded so nobody chases it.** The transfer sends this Windows checkout's
bytes, so 108 files on the box now have CRLF line endings where they previously had LF. Python is
indifferent and the checksum mismatch this produces against an LF-normalised comparison is *not*
content drift. Anything that later diffs the two must normalise line endings first or it will report
108 phantom differences.

### In plain terms

The copy of the code on the GPU machine was not the code in the repository. I had been keeping it up
to date only for the files I personally edited, so every fix made by anyone else was missing — and the
one that mattered was the fix that lets the negative control run at all.

It did not announce itself as "you are running old code". It announced itself as an experiment
failing, and the failure looked like a data problem.

### Meaning for the claim

**Nothing measured on that workspace can be assumed to have used current code**, and that includes the
D1-B training runs, the momentum sweep and the gate variance study. In practice the risk is low: those
depend on `training.py`, `runner.py` and `model.py`, which I edited and therefore kept current. But
"in practice the risk is low" is an argument, not a check, and the correct response is that any number
promoted to the paper should be recomputed against a workspace verified equal to HEAD. The rank
recomputation another agent is running already does exactly this for the rank numbers.

The audit's A3 result, once it lands, is computed on the corrected code.

### Files / commits

- `v2/research/rebase/d2_compare.py` — `_targets`, the group-selection fix
- `~/morpheus-rebase-d1` — now a full mirror of HEAD, committed locally so drift shows as dirty
- Related operational failures of the same family: `operational_shared_box_rules_20260804T0730Z.md`
