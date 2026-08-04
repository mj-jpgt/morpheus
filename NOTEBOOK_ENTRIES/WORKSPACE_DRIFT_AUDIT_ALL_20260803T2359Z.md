## 2026-08-03 23:59 UTC — ALL SEVEN GPU workspaces have drifted from HEAD. Per-workspace table, for the owners to act on

**Logged:** 2026-08-03 23:59 UTC. **How obtained:** `/tmp/audit_workspaces.py` and `/tmp/audit2.py`
on the box, comparing every tracked file against an LF-normalised manifest of HEAD
(`/tmp/head_manifest.txt`, 393 files). **Read-only. No other agent's workspace was modified** — see
"why nothing was fixed for you" below.

### The table

Code = `v2/ tests/ src/ configs/`, `.py|.json|.yaml`. Docs (`NOTEBOOK_ENTRIES/`, `paper/`) excluded:
a missing notebook entry is harmless, a stale `spectral.py` is not.

| workspace | real path | code match | **code DIFFER** | code missing | verdict |
|---|---|---:|---:|---:|---|
| `~/ws` | `morpheus-rebase` | 161 | **15** | 21 | **worst** |
| `~/ws_d1` | `morpheus-rebase-d1` | 197 | **0** | 0 | clean (synced 23:40) |
| `~/ws_d3` | `morpheus-rebase-d3` | 182 | **12** | 3 | drifted |
| `~/ws_k` | *(mine, throwaway)* | 165 | **12** | 20 | disposable |
| `~/ws_p1` | `morpheus-rebase-p1` | 167 | **15** | 15 | drifted |
| `~/ws_rank` | *(rank recompute)* | 165 | **12** | 20 | drifted |
| `~/ws_spatial` | `morpheus-rebase-spatial` | 170 | **15** | 12 | drifted |

### Stale files that compute numbers

| workspace | stale / missing among `calibra/`, `d2_compare`, `spectral`, `training`, `runner`, `paired_bootstrap`, `d1_audit` |
|---|---|
| `~/ws` | DIFFER `calibra/__init__`, `claim_guards`, `e0_basis_transfer`, `run_calibra`, **`spectral`**, `d2_compare`, `runner`, `training`; MISSING `e0_proliferation_stratified`, `gene_label_shuffle_control`, `hest`, `hest_build`, `track1_battery_ledger`, `d1_audit`, `test_d2_compare_deflation` |
| `~/ws_d3` | DIFFER `calibra/__init__`, `e0_basis_transfer`, **`spectral`**, `d1_audit`, `runner`, `training` |
| `~/ws_p1` | DIFFER `calibra/__init__`, `claim_guards`, `e0_basis_transfer`, `run_calibra`, **`spectral`**, `d2_compare`, `runner`, `training`; MISSING 5 |
| `~/ws_rank` | DIFFER `calibra/__init__`, `claim_guards`, `e0_basis_transfer`, `run_calibra`, `d1_audit`, `d2_compare`, `runner`, `training`; **`spectral` is CURRENT** |
| `~/ws_spatial` | DIFFER `calibra/__init__`, `claim_guards`, `e0_basis_transfer`, `run_calibra`, **`spectral`**, `d2_compare`, `runner`, `training`; MISSING 3 |
| `~/ws_d1` | none |

### What this did and did not cost

**`~/ws` is where P2's metric scripts ran** (`p2_rank_variants.py` hardcodes
`sys.path.insert(0,"/home/ubuntu/ws")`), and its `spectral.py` **predates the rank canonicalisation
entirely** — no `CANONICAL`, no `RANK_VARIANTS`. Every P2 rank number was therefore computed by a
different function object than the one now in the repository.

**The numbers nonetheless reproduce exactly.** Recomputed against `~/ws_d1` (verified 0 differ,
0 missing):

| quantity | recorded | recomputed | |
|---|---|---|---|
| effective rank, arm share | 34.5% | 34.5% | ✓ |
| effective rank, F(3,8) | 1.41 | 1.41 | ✓ |
| RankMe, arm share | 29.1% | 29.1% | ✓ |
| ground-truth CCA, arm share | 98.0% | 98.0% | ✓ |
| ground-truth CCA, F(3,8) | 128.20 | 128.20 | ✓ |

They reproduce because the pre-consolidation `effective_rank` is numerically the same function as the
new `CANONICAL` — Roy & Vetterli, order 1, centred, rows at own norms. The consolidation refactored
and added variants without moving the default. `residualise.py` is byte-identical across all
workspaces, so the channel side was never exposed at all.

**That is luck, not design.** The consolidation introduced three named variants; had it chosen a
different default, every P2 rank number would have shifted silently and the reproduction check is the
only thing that would have caught it.

**`~/ws_rank` has a current `spectral.py`**, so the rank recomputation's *definition* was right. But
that workspace has 12 other stale files including `d2_compare.py`, so the assurance covers the rank
definition and nothing else computed there.

### Why nothing was fixed for you

Several of these workspaces have live jobs. Overwriting code under a running process produces a
mid-run failure that presents as a data problem — the exact class this project has spent two days
cataloguing (a harness that was not the runner; a key that returned `[]`; an idle-check that fired in
the gap; a sync that measured the wrong delta). Measuring the extent and leaving the timing to each
owner is the correct division of the work.

**To fix one, do not sync a diff.** That is the bug: `git diff <my-last-commit>..HEAD` is the set
changed since *your* last commit, not the set differing from *your workspace*, so anything fixed
before that range is stale in both and invisible. Ship every tracked file:

```
git ls-files > /tmp/tracked.txt
tar -cf - -T /tmp/tracked.txt | ssh box 'cd ~/<workspace> && tar -xf -'
```

then `git init` + commit inside the workspace so future drift shows as a dirty tree.

### Two workspaces nobody was tracking

`~/ws_k` is **mine** — a throwaway created repeatedly as `rm -rf ~/ws_k && cp -r morpheus-rebase-d1
~/ws_k/morpheus` to test modified files without touching the live one. **It has produced no quoted
number** and can be deleted at any time. Its drift is expected and irrelevant.

`~/ws_rank` is not mine. It appears to be the rank-recomputation workspace. Its `spectral.py` is
current, which is the file that matters for those numbers, but its owner should confirm what else it
produced.

### In plain terms

Every copy of the code on the GPU box is out of date except the one I fixed an hour ago, and the copy
that produced the paper's rank numbers was the most out of date of all — running a version of the
rank function from before it was standardised.

The numbers happen to be identical anyway, because standardising the function did not change what it
computes by default. If it had, nothing would have told us.

### Files / commits

- `/tmp/head_manifest.txt`, `/tmp/audit_workspaces.py`, `/tmp/audit2.py` on the box
- `~/e0_run/p2_*.py` — the P2 metric scripts, **not vendored into the repository**
- Prior: `workspace_drift_stale_code_20260803T2350Z.md`
