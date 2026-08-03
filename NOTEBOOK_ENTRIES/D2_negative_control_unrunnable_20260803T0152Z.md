## 2026-08-03 01:52 UTC — the D2 negative control could not be selected at all (`_targets` bug)

**Logged:** 2026-08-03 01:52 UTC. **How obtained:** reading `v2/research/rebase/d2_compare.py::_targets`
while preparing the stratified readout; verified against `~/e0_run/data/frozen_rna_targets.npz` on
the Lambda box.

### Technical
`_targets` applied the `RANDOM_CONTROL__` name-prefix drop **unconditionally**, then intersected it
with the requested `--target-groups`:

```python
keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
if groups:
    keep &= np.isin(available, list(groups))
```

`--target-groups random_control` therefore selects the intersection of "is not a random control" and
"is a random control" — the empty set — and dies on the next line with
`target-group selection ['random_control'] left no targets`. The negative control that the whole
defensibility argument rests on was unrunnable by construction.

Fix: when `--target-groups` is given, select purely by group; the prefix drop remains the default
"all non-control" behaviour.

```python
if groups:
    keep = np.isin(available, list(groups))
else:
    keep = ~np.char.startswith(names, "RANDOM_CONTROL__")
```

This is provably behaviour-preserving for every non-control selection: in the frozen artifact the
prefix mask and the `random_control` group label are in **exact** correspondence — asserted on the
box, `(prefix == group).all() == True`, 90 and 90. Group census of `frozen_rna_targets.npz`:
180 targets = `hallmark_in_training` 50 + `heldout_pathway` 24 + `immune_tme` 8 + `tumour_state` 8 +
`random_control` 90. So the unrestricted 90 and the untrained 40 are unchanged by this edit.

### In plain terms
The switch that was supposed to let us score the models on deliberately meaningless targets — the
check that tells us whether the measuring instrument invents signal out of noise — was wired so that
asking for those targets asked for "the targets that are both random and not random". It always
returned nothing and crashed. Nobody would have noticed except by trying to run it.

### Meaning for the claim
Directly load-bearing. Without this the paper's "we ran a negative control" sentence could not have
been written truthfully. It also means **no previous D2 result was ever checked against a negative
control**, so the −0.10 headline has never had that guard applied to it.

### Files / commits
- `v2/research/rebase/d2_compare.py` — local `d4aaf96`, Lambda `3853d40`
