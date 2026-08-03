# The claim guards were documentation, not guards — `validate_claim` read evidence from nowhere. They now read a provenanced, checksummed evidence file, and the first thing the new check did was refuse a discharge whose notebook entry was missing.

**Logged:** 2026-08-03 15:45 UTC. Requirements pre-registered in
`NOTEBOOK_ENTRIES/e0_proliferation_stratification_preregistration_20260803T1450Z.md` §5, committed
`5f331f6`, **before** the code was written.

**How obtained:** code + tests, CPU only, thread-capped. Full suite **302 passed, 0 failed**.

---

### Technical

#### The defect

`claim_guards.validate_claim(claim: dict)` is pure logic over a dict somebody hands it. Grepping the
repo for who builds that dict returns **nothing in production** — only the tests. The project's record
of whether E0 was publishable lived at `tests/test_claim_guards.py:135` as:

```python
e0 = {"kind": "transfer", "proliferation_controlled": False,
      "platforms": ["perturb_seq_crispri", "perturb_seq_crispri"]}
```

So the guard could not *catch* anything. It restated a verdict a human had already typed, and
"discharging a blocker" meant editing a test to say the blocker was discharged. That is the same
failure mode this project keeps finding in its experiments — **a check that looks like it ran and
didn't** — sitting in the machinery built to prevent exactly that.

#### What was built

`v2/research/rebase/nature/claim_evidence.json`, read by
`claim_guards.validate_recorded_claim(name, path, repo_root)`. Each evidence field is a record, not a
bare value:

```json
"proliferation_controlled": {
  "value": true,
  "run":    "~/e0_run/d3/e0_prolif/e0_proliferation_stratified.json",
  "entry":  "NOTEBOOK_ENTRIES/e0_proliferation_stratification_result_20260803T1530Z.md",
  "commit": "bfe4569",
  "how":    "...", "scope": "..."
}
```

Against the six pre-declared requirements:

| # | requirement | how |
|---|---|---|
| 1 | checksummed evidence file, provenance triple per field | `run` + `entry` + `commit`, plus `how`/`scope` prose |
| 2 | missing/unreadable/mismatched ⇒ **inadmissible**, never permissive | `load_claim_evidence` returns `{}` + a note; `validate_recorded_claim` returns every blocker for the kind |
| 3 | **evidence with no provenance is not evidence** | `_resolved_value` requires a non-empty `run`, an `entry` that **exists on disk**, and a `commit` matching `^[0-9a-f]{7,40}$`; otherwise the field is reported ABSENT and its blocker fires |
| 4 | fixture test rewritten to assert **logic** | the pinned-verdict test is gone; 25 tests now exercise each blocker firing and clearing, plus every failure mode above |
| 5 | real claim state moves out of the test | two records: `E0_basis_transfer_K562` and `morphology_to_pbs_axis_legibility` |
| 6 | a test that the guard **still bites** | `test_the_recorded_e0_claim_is_still_inadmissible` |

The `sha256` is over the canonicalised evidence block excluding itself. **It is tamper-evidence, not
tamper-proofing, and is documented as such in the file** — anyone can recompute it. What it buys is
that a value cannot change *silently*: the digest must change too, which is visible in a diff.

#### The refactor's own falsifier, pre-declared and now asserted

> *"If after the change a claim can be made admissible by editing only the evidence file, with no
> analysis run and no provenance, the refactor has failed and must be reverted."*

`test_REFACTOR_FALSIFIER_a_claim_cannot_be_discharged_by_editing_values_alone` writes every discharge
field as `true` with an **honestly recomputed digest** and no provenance, and asserts the claim is
still inadmissible with `{proliferation_deflation, single_platform}` intact. It passes.

#### The new check earned its keep immediately

On its first run against the real file, `test_the_shipped_evidence_file_is_internally_consistent`
failed with:

```
E0_basis_transfer_K562: proliferation_controlled: notebook entry
'NOTEBOOK_ENTRIES/e0_...T1530Z.md' does not exist; treated as absent
```

The workspace copy had no `NOTEBOOK_ENTRIES/`, so the discharge was refused and
`proliferation_deflation` was reported as still blocking. The old system could not have noticed —
there was nothing to notice with. Under the old fixture, the same evidence would have been a `True`
in a test file with no way to tell whether the entry behind it existed.

#### Current recorded state

| claim | kind | admissible | blockers |
|---|---|---|---|
| `E0_basis_transfer_K562` | transfer | **No** | `single_platform` |
| `morphology_to_pbs_axis_legibility` | legible_axis | **No** | `composition_attribution`, `no_external_cohort` |

`proliferation_deflation` is discharged for E0 by
`e0_proliferation_stratification_result_20260803T1530Z.md`; `purity_confound` is discharged for the
legibility claim by `d3_purity_result_20260803T1330Z.md`. **Both claims remain inadmissible.** The
refactor did not loosen anything.

### In plain terms

We had a system that was supposed to stop us publishing a claim before the right controls existed. It
turned out it never looked at anything. It was a list of rules plus a note, written by hand, saying
which rules had been satisfied — and the note lived inside a test. To mark a caveat as handled you
edited the test that asserted it wasn't.

Now the rules read a file that says, for each control, which run produced it, which write-up describes
it, and which commit it landed in. If any of those three is missing or points at something that isn't
there, the control counts as not done. So you can no longer discharge a caveat by typing the answer
you want; you have to have run something and written it down.

The first time we ran it, it caught us: a write-up referenced in the file wasn't present in that
working copy, and it refused the discharge.

### Meaning for the claim

- Nothing became admissible. Two claims are recorded and both are still blocked.
- `proliferation_deflation` moves from "undischarged" to "discharged for E0, with provenance and a
  declared scope" — recorded as **data** rather than as an edit to a test assertion.
- **Honest limits of the new system**, so it is not oversold in turn:
  1. It checks that provenance *resolves*, not that it is *true*. A field can name a real run, a real
     entry and a real commit whose content does not support the value. Only review catches that.
  2. `commit` is format-checked, not verified against git history — the guard must run without a git
     worktree.
  3. Nothing yet *writes* the evidence file from an analysis. It is still authored by hand; the
     analyses do not emit their own records. That is the next step and it is not done.
  4. Only two claims are recorded. Every other claim in the project remains unrecorded, and an
     unrecorded claim is inadmissible by default — correct, but it means coverage is thin.

### Files / commits

- `v2/research/rebase/nature/claim_evidence.json` (new)
- `v2/calibra/claim_guards.py` — `CLAIM_EVIDENCE_PATH`, `evidence_digest`, `load_claim_evidence`,
  `validate_recorded_claim`
- `tests/test_claim_guards.py` — pinned-verdict test replaced by 8 evidence-layer tests
