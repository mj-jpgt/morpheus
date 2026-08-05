# PROJECT GUIDE — read this first, every time

**This document supersedes the entire `HANDOFF_*.md` chain as the entry point.** Those files (`HANDOFF_BUILD_AGENT.md`, `HANDOFF_EXPERIMENTS_NOW.md`, `HANDOFF_GATES.md`, `HANDOFF_PHASE_D.md`) cover the E-series and Phase D build history, which is now folded into the paper drafts below. They are kept for archaeology and are individually still correct about what they describe — do not treat them as the active task list. **This file is not another link in that chain. It does not get superseded by a newer handoff file — it gets edited in place, forever, by whoever's work changes what it says.** See §5.

If you are an agent dropped into this repository with no other context, read this file, then `NOTEBOOK.md`'s most recent ~30 entries, then the paper draft(s) relevant to your task. That is enough to be useful without breaking something.

---

## 1. What this project is

MORPHEUS/biorag is a multimodal tumour representation project: whole-slide H&E images, bulk RNA, and (in progress) copy-number, mutation, and protein data, trained with contrastive self-supervised learning (paired InfoNCE with a MoCo-style momentum key encoder — see `v2/losses.py`, `v2/training.py`) into a shared representation, held to a confound-adjustment measurement instrument called **CALIBRA** (`v2/calibra/`) that asks, rigorously, whether the representation carries molecular information a tumour's cancer type and originating site cannot already explain.

Four papers, deliberately kept separate (do not merge P1/P2, do not split gates into their own paper — both were explicit decisions):

- **P1 — CALIBRA** (`paper/P1_CALIBRA_DRAFT.md`): the instrument itself, and its headline result — a morphology→molecular channel that survives the strongest confound adjustment that is mathematically possible, validated out-of-sample, replicated in an external cohort (ALCHEMIST). The project's strongest, most-tested claim.
- **P2 — effective rank** (`paper/P2_RANK_DRAFT.md`): a negative methods result. Effective rank's usability as a label-free selection signal is conditional on the co-trained view and the statistic used to read it — and, newly established, on the reference standard it's meant to substitute for having the same instability, not a defect unique to rank.
- **P3 — perturbation-basis / causal attribution**: an interventional gene-perturbation dictionary lost, cleanly and repeatedly, as a *competing basis* against plain PCA. Pivoted to a working, different contribution — using the same resource to attach certified causal names to axes a representation already has.
- **P4 — the promptable multimodal system**: certifies before answering. Currently answers close to nothing, for a mix of a real fixable engineering gap (no inductive out-of-sample adjustment existed until today) and a real ceiling (mostly only answers what it was explicitly trained on) that is under active, separate investigation.
- **P5 — a discovery engine** (`paper/P5_DISCOVERY_PLAN.md`, planning stage, no measurement yet): sweeps the representation against a large candidate space of gene/pathway/mutation targets to surface novel, statistically certified, replicated, mechanistically-grounded hypotheses with a proposed validation experiment attached — a generation system, not a retrieval one. Kept distinct from P4 on the same separation logic as P1/P2. The hard, largely unbuilt part is filtering out rediscovery of already-known biology; see the plan's §1.

Branch: `research/rebase-vision`. Do not start work from `main` or `v2` — they predate all of the above.

---

## 2. How to work on this project — the methodology, and why each rule exists

These are not style preferences. Every one of them exists because its absence already produced a real, documented failure on this project. Cite the incident if you want to know why a rule is here; they're in `NOTEBOOK_ENTRIES/`.

1. **Predeclare before you measure.** Write what result would mean what — including what would make you *distrust* a favourable result — to `NOTEBOOK_ENTRIES/PREDECLARED_<name>_<UTC-timestamp>.md`, committed **before** the measuring code runs, not after. *Why:* a "centring explains the dissociation" theory once fit its own math exactly (predicted 2.106 against 2.107 measured) and was still wrong, because the motivating facts were a coincidence of which blocks were compared. Predeclaring the falsifier in advance is the only thing that catches that kind of self-deception.
2. **Report bad news first.** Every notebook entry, every agent report, every summary to the user leads with what's wrong, what failed, or what's still unresolved — not with what succeeded. *Why:* a headline number (`0.463 → 0.035`) sat quoted in the paper for weeks because nobody's first instinct was to ask what could be wrong with a favourable-looking result.
3. **Push every favourable result until it breaks, or until you're confident it can't.** More repeats, a harder statistic, more partitions, an out-of-sample test, a bigger discovery fold. Do not accept n=5 same-seed repeats as a floor if n=10 is affordable. *Why, concretely:* the P2 momentum claim (`m=0.999` over `m=0.99`) cleared its floor at n=5 and **failed** at n=10 — caught only because someone was told to push on it rather than bank the pass.
4. **A negative result or a fair-test loss is a checkpoint to redesign, not a conclusion to stop at.** Report it plainly, log it, and then propose the next construction. *Why:* PBS losing to PCA didn't end P3 — it produced the causal-attribution pivot, which is a better contribution than the one that failed.
5. **Reuse canonical functions. Never write an inline reimplementation of a statistic that already exists in `calibra`.** *Why:* this project has had five separate silent statistic-substitution bugs — an inline formula standing in for an import, computing something subtly different under the same name. An AST-scan test (`v2/tests/test_effective_rank_canonical.py` and siblings) now fails any new inline `linalg.svd`/rank computation that doesn't import the canonical function; do not work around it, import instead.
6. **Every claim ships with its counterfactual.** A positive control that must clearly pass, a must-fail control that must clearly fail, a matched random-direction or permutation null, and — for any certificate or gate — a planted-signal ladder proving the instrument actually discriminates rather than rubber-stamping everything. *Why:* P4's confound certificate was validated this way (planted site/cancer/noise codes at four strengths) specifically because a certificate that always passes is worse than no certificate.
7. **Verify a claim traces to a real artifact with a real hash before trusting it.** *Why:* three separate incidents of the same filename pointing to different files with different answers (`d2_h_seed42.npz`, `diagnostic_full_seed42.npz`), and one headline number that, on audit, traced to no artifact and no code path at all.
8. **Every agent logs everything to the notebook, every time, no exceptions.** Write a full entry to `NOTEBOOK_ENTRIES/<descriptive-name>_<UTC-timestamp>.md` for anything you predeclare, measure, build, or decide — including negative results, things you tried and abandoned, and bugs you found in someone else's work. **Never edit `NOTEBOOK.md` directly** — the coordinating session merges the index from entries. A result that isn't in the notebook did not happen, for the purposes of this project.
9. **Corrections are appended, not silently rewritten.** History stays visible. See the `CORRECTION APPENDED` block convention (e.g. `t13_adjusted_certificate_and_p6`) and the numbered-status-item-append convention in `P2_RANK_DRAFT.md`'s status block — never renumber or delete a prior item, add the next one.
10. **Shared-tree git discipline.** Several agents edit this repository concurrently, often the same files. `git add` only the specific paths you changed — never `git add -A` or `git add .`, which has swept another agent's in-progress work into a commit before. `git pull --rebase` before pushing; if that fails because another agent has **uncommitted local changes** (not a rebase conflict), do not touch or stage their files — run `git fetch && git status -sb`, and if you're only "ahead" of origin (not diverged), a plain `git push` works without needing a clean tree.
11. **Workspace sync uses a manifest, never a git-diff-based sync.** `git diff <my-last-commit>..HEAD` measures the wrong delta and has caused real drift across seven box workspaces at once. Use `--workspace-manifest`, kept outside the workspace it verifies.
12. **`git archive` honours `core.autocrlf`** — a file-count check after syncing is not sufficient; verify by content hash (blob SHA or md5 manifest), not count.
13. **Watch for self-matching process patterns.** `pgrep -f`/`pkill -f` can match your own polling command's argv (e.g. `pgrep -f "pytest morpheus"` matching a `bash -c "... pytest morpheus ..."` wrapper around itself). Prefer a PID file or a more specific pattern.
14. **Cap BLAS threads on the shared box.** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` — uncapped, one 6-second test took over 25 minutes from thread oversubscription with other agents' jobs running concurrently.
15. **Check GPU occupancy before training; queue, don't force.** `nvidia-smi` compute-app count and memory, not a process-name guess. Launching into a saturated card risks OOM-ing another agent's multi-hour run. If queued, say so and report what's completable without the GPU in the meantime.
16. **A per-axis or per-arm-only measurement is not the whole story.** Measuring only the stable arm, only one partition, only one reading step, or only one statistic has repeatedly produced a flattered result on this project (the P2 headline floor measured on the stable arm only; P1's inductive channel first reported on one partition; audit verdicts that flip four times across six reading steps). Default to sweeping the axis you'd otherwise fix, or say explicitly why you didn't.

---

## 3. Current state and remaining build plan, per paper

*(Snapshot as of this section's last edit — see the top of `NOTEBOOK.md` for anything more recent than this document.)*

### P1 — CALIBRA (strongest, most tested)
**Solid:** confound adjustment verified on two confounders and provably bounded (saturated cell design upper-bounds any conditional-mean method); channel survives out-of-sample adjustment across 12 independent discovery/exposure partitions, not just one; replicates externally on ALCHEMIST (R=1.110, n-matched); a general methods finding along the way (cross-fitted adjustment manufactures apparent confound recovery from pure noise) that reaches beyond this project.

**Open:**
- ~~Citation verification pass~~ **DONE 2026-08-05** (`NOTEBOOK_ENTRIES/p1_completeness_pass_two_wrong_citations_and_the_external_cohort_denial_20260805T0430Z.md`). All eleven `[UNVERIFIED]`/`[CITATION NEEDED]` markers resolved against live records. **Two citations were wrong**, both author-name errors on real papers: the random-signature control was "Venet, **Dhanasekaran & Sotiriou**" (actually Venet, **Dumont & Detours**) and ComBat was "Johnson, Li & **Rabinowitz**" (actually **Rabinovic**). Both corrected in the draft. Residual: the *page* in Muirhead/Anderson for the N−R result — needs the books, not an API.
- ~~Full read-through for internal consistency~~ **DONE 2026-08-05**, same entry. Six defects fixed in `P1_CALIBRA_DRAFT.md` and four in `P1_FIGURES.md`. The one that mattered: the abstract, §1.3 and §5.1 all asserted "**No external cohort has been through the instrument**" — false since ALCHEMIST landed 2026-08-04 21:15. Now scoped to what is actually true (the *floors* have never been measured outside TCGA; the *channel* has been replicated).
- The dangling `morphology_to_pbs_axis_legibility` claim in `claim_evidence.json` is moot now that PBS lost to PCA — needs an explicit withdrawal/supersession note, not silent deletion. **Still open. Verified 2026-08-05 that no paper draft references this claim record by name**, and that P1 §4.13 already carries the table that moots it. Safe mechanism, checked against the guards: a **claim-level** `"status": "withdrawn"` / `"superseded_by"` key (sibling of `kind`/`description`/`evidence`) is ignored by `load_claim_evidence` and is **not** covered by `evidence_digest`, so it will not trip `test_the_shipped_evidence_file_is_internally_consistent`. Do not put it inside `evidence`.
- `no_external_cohort` deliberately stays undischarged for any **per-axis** claim (it gates `legible_axis`/`gene_attribution` kinds specifically) — ALCHEMIST replicated the aggregate channel, not any named axis. Do not discharge this blocker with the ALCHEMIST result without a per-axis external replication. *P1 §1.3 and §5.1 now say this in the paper, so the reason is no longer only in a notebook entry.*
- **New, needs measurement not writing:** `TRACK1_NEGATIVE_CONTROLS.md` §T1.3 states the adjusted site-LDA drop as "21–45×"; recomputed from the six rows of its own table the range is **26–45×** (no row gives 21). Corrected in P1; the source result file still says 21–45×.
- **New:** F2 panel (f) (the twelve-partition inductive retention strip) is specified in `P1_FIGURES.md` but marked `NEEDS EXTRACTION` — the 24 retention values are only in a notebook-entry table, not in a plot-ready file.

### P2 — effective rank (negative result, actively being hardened)
**Solid:** the view-conditional instability is real and now attributed correctly — to the view, not the metric (a labelled linear probe, the actual reference standard, shows the identical instability). A 62-row self-audit with per-statistic, per-block, per-step floors. Citations verified (13/13 arXiv IDs real).

**Open:**
- The paper's headline floor (3.295×, canonical R1 on the exported block) was measured on the **stable** training arm only — the paper's own analysis says that flatters a floor by ~2×. An unstable-arm measurement is training now; **do not treat 3.295× as settled until that result lands.**
- §5.2's `m=0.9` dip (2.23 against `m=0`'s 2.81, one seed each) is unresolved and sits in the one interval the momentum grid doesn't cover — flagged, not measured.
- §4.1a's floor table needs a generator change (`p2_floor_audit.py` extended to take `P2_LABELLED_PROBE.json`) before the probe rows can join the main rendered table instead of sitting in a separate one.

### P3 — perturbation-basis / causal attribution (pivoted, live contribution identified)
**Solid:** PBS-as-competing-basis is dead across five honest attempts (original, joint CCA, cross-line consensus, domain-adapted) — stop pursuing that framing. The pivot — attributing already-legible axes to causal gene-perturbation signatures — has a real four-condition certificate; 29/128 PCA axes pass, landing on recognizable complexes (ribosome biogenesis, mitochondrial nucleoid, cohesin).

**Open (agent running as of this writing):**
- Whether a basis chosen to be attributable (rather than variance-maximizing) certifies substantially more than 29/128 — tests whether 29/128 is a ceiling on the biology or an artifact of using raw PCA axes.
- Whether the 29 certified axes disproportionately carry the actual morphology→molecular channel, or whether certification is (as suspected) roughly orthogonal to which axes matter for prediction.

### P4 — promptable system (earliest stage, most open)
**Solid:** the five-condition certificate discriminates correctly (validated on a planted-axis strength ladder). The inductive (out-of-sample) adjustment operator now exists and is wired in.

**Answered 2026-08-05, and it's the negative branch:** both generalization tests landed. **(A) — no signal, not a readout-shape problem.** Held-out pathways are answered by the incumbent rule *less often than a matched random gene set* (8.3% vs 23.3%); a composed multivariate readout raises both together (29.2% vs 36.7%) because the rise is a capacity artifact, not specificity — nine grading rules, nine failures to beat the matched control. **The causal-name bridge (P3→P4) also fails**: certified names carry no predictive power for what an axis reads (Spearman +0.013, p=0.334), and the mechanism is now understood — the ten strongest reads in the whole representation are all immune targets, none of the 29 certified axes is named for anything immune, because the K562/RPE1 causal dictionary can only ever name cell-autonomous processes, not tissue-composition ones. Full detail and the P5 implications: `paper/P5_DISCOVERY_PLAN.md` §1.1/§2.1 — this also surfaced a real gap in CALIBRA's own floor/null machinery (doesn't control for readout capacity) that belongs in P1 too, not just here.

**Open:**
- Multi-modal fusion architecture is **unstarted**. Current training is strictly two-tower (WSI↔RNA) InfoNCE; five modalities with heterogeneous missingness (SNV ~95% coverage, CNV ~96%, RPPA only ~71%) don't extend cleanly to pairwise contrastive. Candidate direction: shared latent with masked-modality training (missingness at inference = masking at training), to be decided against measurement (held-out-modality reconstruction, whether the channel survives) rather than architectural preference — not yet built or scoped in code.
- SNV modality is built and shows structure after adjustment, but not yet distinguished from the lineage/cancer-type residual found elsewhere on this project — needs its own channel test before it's a claim.
- CNV and RPPA proteomics are scoped (data located, coverage measured) but not built.

---

## 4. Where things live

- `v2/calibra/` — the instrument: adjustment (`residualise.py`, `inductive_adjustment.py`), measurement (`spectral.py` — the **only** place rank/CCA statistics are defined, everything else imports from here), certification (`confound_certificate.py`, `claim_guards.py`), and honest-null machinery (`calibration.py`, `honest_metrics.py`).
- `v2/losses.py`, `v2/training.py` — the SSL training objective (paired InfoNCE + momentum queue).
- `NOTEBOOK.md` — the merged, indexed log. Read the most recent entries for anything more current than this document.
- `NOTEBOOK_ENTRIES/` — where every agent writes; the primary source of truth for what actually happened and why.
- `paper/` — the four drafts, plus `*_FIGURES.md` plans.
- `v2/research/rebase/nature/claim_evidence.json` — the machine-checked claim admissibility registry.
- `v2/tests/` — includes several guard tests that exist specifically to catch the failure modes in §2 (statistic duplication, workspace drift, artifact provenance). Do not weaken or work around one of these; if it's blocking you, the block is very likely correct and the fix is elsewhere.

**Two test-infrastructure hazards, both real incidents:**
- `pytest v2/tests tests` from the repo root fails with ~70 collection errors — the repo must be importable **as the package `morpheus`**, which the checked-out directory name (`morpheus-rebase`) is not. Run from wherever your workspace convention puts `morpheus` on `PYTHONPATH`, not naively from the checkout root.
- **A copied (not symlinked/junctioned) workspace silently reports a stale suite as current.** One agent's `scratchpad/ws/morpheus` was a stale *copy* rather than a link — 47 test files instead of 58 — and its "517 passed" was caught only because it happened to exactly match a count from the *previous day*. If a suite count looks suspiciously familiar, verify the workspace is actually current (file count, or better, blob-hash manifest) before trusting the number.

---

## 5. Keeping this document current — not optional

**Any agent whose work closes an item in §3, opens a new one, or changes the build plan must update the relevant part of this document before finishing** — in addition to, not instead of, its `NOTEBOOK_ENTRIES/` entry. The notebook is the detailed record of what happened; this document is the current summary someone can read cold. They serve different purposes and both are required.

If your change is small (a line in §3's open-items list), edit it directly. If it's structural (a new paper direction, a methodology rule earned from a new incident), add it in the matching style — cite the concrete incident that motivated it, the way every rule in §2 does. Do not let this document drift back into being a snapshot of one day; it is meant to always be readable as "the state of the project right now."
