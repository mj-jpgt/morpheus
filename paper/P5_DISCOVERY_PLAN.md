# P5 — a discovery engine, not a query system: generating novel, testable biological hypotheses at scale

**Compiled:** 2026-08-05 · **Branch:** `research/rebase-vision` · **Status:** planning document, no
measurement yet at pipeline scale (see §2 for one real result that already lands on this plan
directly). Written in response to a direct instruction: build a plan for genuine in-silico
biological discovery — novel, cutting-edge hypotheses the field doesn't already have an answer for,
generated at mass scale from the existing multimodal representation, ranked, and handed off with a
concrete path to in vitro/in vivo validation. Not a re-derivation of known biology; not a per-patient
answer to a pre-set question.

**One caveat stated up front rather than buried:** the specific external methodology referenced as a
scale exemplar could not be verified — the session's web search budget was exhausted when this was
written. Nothing below depends on that reference; it's built from the requirement as stated
(mass-scale generation, then validation), and can be revised once that comparison is actually checked.

---

## 0. Why this is P5, not P4 — recommendation, not a unilateral decision

P4 is a **retrieval system**: certify whether a *given* answer to a *given* query about a *given*
patient is admissible. Its unit of output is "yes, with a certificate" or "no, and here's why not."

What's being asked for here is a **generation system**: sweep the representation against a large,
open candidate space of gene/pathway/mutation targets across the cohort, and surface the *small
number* of statistically real, novel, mechanistically plausible associations nobody has looked for
yet. Its unit of output is a ranked shortlist of falsifiable hypotheses with a proposed experiment
attached. The audience is different (a wet lab or a collaborating biologist, not a clinician at the
point of care), the validation standard is different (an experiment outside this project entirely,
not a held-out patient), and the failure mode is different (a false positive here wastes someone's
bench time, not a clinical decision).

That's a difference in kind, not just scope, and this project's own precedent (P1/P2 kept separate
despite sharing an instrument; gates never got their own paper) favours keeping this as its own
paper. **Recommendation: P5, sharing infrastructure with P1/P3/P4 but reported as a distinct
contribution.** Overridable — say so if you'd rather fold it into P4 as a section.

---

## 1. What kind of question this system is actually suited for — and what it isn't

This section exists because "find novel biology" is not itself a scope — it's crowded, and most of
it is not where this system's advantage is. Stated plainly, then argued:

**Not well-suited:** single-gene mutation-from-image discovery (a crowded field already — Kather,
CHIEF, Prov-GigaPath, all cited in P1 — and one this project's own audits show is often secretly
cancer-type detection, not the claimed biology). **Not well-suited:** molecular-only discovery with
no image involved (pure genomics pipelines like DepMap have no confounding from tissue to strip and
are simply better positioned there).

**What's actually ours is the intersection of four things, and each maps to a piece of
infrastructure that already exists:**

1. **A confound-adjustment standard stricter than the field's own published work supplies.** Most
   H&E→molecular papers don't strip cancer-type and site the way this project does — provably bounded,
   tested inductively, validated across 12 partitions (P1). That licenses a distinctive question:
   *of the morphology-molecular links already claimed in the literature, which survive real
   deconfounding, and which are secretly lineage detection?* — separating real from spurious among
   *existing* claims, with an instrument stricter than what produced them.
2. **The channel is multivariate, not per-gene** — P1 already found 76–82% of the per-target channel
   is reproduced by covariate-matched random gene sets. That's not a defect to route around; it's a
   direct finding about the honest unit of discovery here: not "does the image show mutation X," but
   *which coordinated multi-gene programs are morphologically visible.*
3. **The causal-attribution bridge (P3) is the genuinely unique asset** — a validated morphology
   channel wired to a Perturb-seq-derived causal dictionary, turning "image correlates with expression
   pattern X" into "perturbing gene G in a standard cell line should reproduce or abolish this visible
   phenotype." Most pathology-genomics work has no causal dictionary in the same coordinate space and
   cannot propose an experiment this specific.
4. **Held-out-whole-cancer-type design points at pan-cancer, lineage-independent programs** —
   because cancer type is explicitly stripped as a confound and generalization is tested across held
   -out cancers, the class of finding this system can honestly reach is *tumour-intrinsic processes
   visible across histologically distinct cancers*, not a marker specific to one cancer's baseline
   composition.

**Put together:** the honestly-scoped target is *coordinated molecular programs that leave a
genuinely non-confounded, cross-cancer morphological signature — and for the subset traceable to
specific gene perturbations, the minimal experiment that should reproduce or erase that visible
phenotype.*

### 1.1 — Two structural limits on that scope, confirmed empirically, not hypothetical

**The dominant morphological signal is very likely tissue composition (immune infiltration,
stroma, necrosis, cellularity) — old, well-known biology, not novel.** `claim_guards.composition_attribution`
has sat undischarged in this project specifically because nobody has shown the channel isn't just
detecting cell composition. This is now more than a theoretical risk: `NOTEBOOK_ENTRIES/` (2026-08-05,
the P4 causal-name-bridge test) found the ten strongest reads in the whole representation are **all
immune targets**, and none of the 29 causally-certified axes is named for anything immune — the best
immune-reading axis (PCA_013) is named for the mitochondrial nucleoid. **Any candidate whose signal
plausibly comes from composition must be filtered out before it's called novel**, and the composition-
attribution control belongs in the funnel's stage 5/6 (§2), not left as an afterthought.

**The K562/RPE1 causal dictionary can only ever ground cell-autonomous processes, not
tissue-composition ones — a structural mismatch, not a data-quantity problem.** Single-gene CRISPRi in
two cultured cell lines cannot produce a signature for a multi-cell-type tissue phenomenon like immune
infiltration; there's no perturbation in that resource that corresponds to "more T cells in the
sample." The same 2026-08-05 test found certified-axis causal names carry **no predictive power** for
what an axis actually reads (Spearman +0.013, p = 0.334; a legibility-matched uncertified control
reads nominally *higher*). Consequence for this plan: **stage 6's mechanistic grounding (§2) is only
available for the cell-autonomous slice of candidates** — the 29 (or however many a better basis
certifies) axes whose domain is things like ribosome biogenesis, mitochondrial function, or cohesin.
Immune/stromal/composition candidates may still be statistically real and worth reporting, but they
cannot honestly carry a K562/RPE1-grounded experiment proposal; a different causal resource (an immune
perturbation atlas, a cytokine-stimulation expression panel) would be needed, and none is currently on
this project.

**Sharpened scope, accounting for both limits:** the strongest, most defensible P5 output is a
**cell-autonomous, non-composition, coordinated molecular program**, novel against the tiers in §3,
mechanistically grounded via the K562/RPE1 bridge, morphologically visible across held-out cancers.
Composition-driven and immune findings can still be reported, but flagged as *not causally groundable
by this project's current resources* and *not novel* by default until proven otherwise — the opposite
of the presumption for cell-autonomous ones.

---

## 2. The discovery funnel — mirrors a wet-lab screening funnel on purpose

Mass scale in, a short, defensible list out. Each stage is a real filter, and — matching this
project's ledger discipline — **every candidate that enters the funnel is tracked through it**,
including the ones that die at each stage, not just the survivors.

| stage | question | cost | reuses |
|---|---|---|---|
| **0 — enumerate** | what's the full candidate space? | cheap, combinatorial | — |
| **1 — coarse pre-filter** | cheap correlation scan, uncorrected, purely to cut compute | cheap | `calibra.spectral` |
| **2 — certify** | confound-adjusted, permutation null, injection-certified floor, **inductive** (out-of-sample, not transductive — see `PROJECT_GUIDE.md` §2 rule 3) | the existing instrument, at scale, **with the capacity caveat in §2.1** | `calibra.*` unchanged |
| **3 — correct for multiplicity** | BH-FDR across the *whole tested space*, predeclared before running | new — nothing on this project has swept an open space before; §4/§90-target work tested pre-chosen targets, not this | — |
| **4 — replicate** | does it hold in an independent cohort (ALCHEMIST; held-out TCGA cancer types)? | the existing infra | `alchemist_*`, `leave_sites_out.py` |
| **5 — novelty filter** | is this already known, at the three tiers in §3? **and** is it plausibly composition-driven (§1.1)? | new, largely unbuilt | — |
| **6 — mechanistic ground** | does an independent functional-genomics resource predict the same direction? does a P3-certified causal name attach? **only available for the cell-autonomous slice (§1.1)** | this project's specific asset | K562/RPE1 Perturb-seq loaders, `causal_attribution.CERTIFICATE` |
| **7 — score and propose** | rank by effect size over floor, replication strength, mechanistic consistency, novelty confidence, and **experimental feasibility**; auto-draft a minimal falsifiable experiment with a stated quantitative prediction | new | — |

**Stage 3 is the statistical crux and deserves its own paragraph.** Every certified result on this
project so far has tested a small, pre-chosen, frozen target list (~90 targets) — the multiple-testing
burden is bounded and known. An open sweep across genes/pathways/cancer-strata is a different
regime — potentially thousands to low-millions of candidate cells before the coarse pre-filter. The
search space size must be **fixed and stated before stage 2 runs**, exactly like every other
predeclaration on this project, so the FDR correction is honest rather than computed against
whatever survived. Two-stage design (FDR-controlled discovery + independent-cohort replication,
mirroring how genomics already does this) is the right shape — a candidate needs to clear *both*,
not either.

### 2.1 — A real defect in stage 2's own instrument, found while testing P4, and it belongs here

The same 2026-08-05 test (`NOTEBOOK_ENTRIES/`, the composed-readout generalization test) found
**CALIBRA's floor and permutation null do not fully control for readout capacity.** A multivariate
readout over more axes inflates *both* the real-target hit rate and the matched-random-control hit
rate together — the floor rose on only 43.3% of targets when capacity increased, and the null actually
*fell* while the statistic rose 2.1×. The only control in that test that correctly absorbed the
capacity increase was a **matched random gene set run through the identical readout**, not the
existing floor or permutation null. **Consequence for this plan: stage 2 must score every candidate
against a matched-random-target control at the same capacity, not only against the permutation floor**
— the floor alone is not sufficient once the readout composes across axes, which any pathway-level
(not single-axis) discovery candidate will. This is a general instrument finding that also belongs in
P1's floor/detection-limit section, not only here.

---

## 3. The actual hard problem: "novel" is the part nothing on this project does yet

Everything else here — certified association, replication, mechanistic grounding — already has
working infrastructure (§4). **Filtering out rediscovery of known biology does not**, and it's the
part that decides whether this is a discovery engine or an expensive way to rediscover MSigDB.

Three tiers of "already known," each needing a different check, from cheap to expensive:

1. **Already in the curated gene-set/pathway resources already on this project** (Hallmark, KEGG
   Medicus, the frozen target tables). Free, already have the data, should be the first filter — if a
   candidate association is just "gene G is a member of pathway P" and P is already a supervised
   target, it isn't a discovery, it's the training signal reappearing.
2. **Already documented in a freely accessible curated database** — Open Targets, DepMap's public
   summary statistics, COSMIC's free tier, GTEx associations. Scriptable, free, not yet integrated
   into any pipeline here.
3. **Already published, not necessarily in a structured database.** PubMed's E-utilities API is free
   and scriptable — a targeted query ("<gene> AND <cancer type> AND (histology OR morphology OR
   imaging)") can surface an existing paper even when no database entry exists. **This tier cannot be
   made airtight.** A negative literature-search result means "not obviously published," not "not
   published" — absence of an indexed hit is evidence, not proof. State this limitation on every
   surfaced hypothesis, not just once in this plan.

None of tier 2 or 3 has been scoped for actual API access, rate limits, or reliability from wherever
this runs (local box vs. the Lambda instance). **That scoping is the first concrete task**, before
any claim is made that the pipeline can tell novel from known — see §6.

---

## 4. What already exists and plugs in directly — this is not starting from zero

- **The instrument.** `v2/calibra/` — confound adjustment (transductive and, as of today, inductive),
  the confound certificate, injection-certified detection/transmission floors, permutation nulls. Zero
  new statistics needed for stages 1, 2, or 4 — but see §2.1's capacity caveat before trusting stage 2
  alone.
- **A second cohort for replication, already wired.** ALCHEMIST (1,106 paired NSCLC patients) — the
  channel already replicates there in aggregate (R=1.110). Per-candidate replication reuses the same
  machinery.
- **A free functional-genomics cross-check, scoped now to its real domain (§1.1).** The K562/RPE1
  Perturb-seq resource (`e0_basis_transfer.py`, `pbs.py` loaders) already on disk. For any
  cell-autonomous candidate gene it covers, stage 6 is a real, in-silico, no-new-data-needed functional
  plausibility check — but it cannot ground tissue-composition candidates; do not force that mapping.
- **A mechanistic-naming certificate, with a measured false-positive risk if misused.** P3's
  `causal_attribution.CERTIFICATE` (four conditions, 29/128 axes certified as of the last basis). A
  certified name does **not** currently predict what an axis reads (§1.1) — use it to add a causal
  story to a candidate the statistics already support, never as a discovery signal by itself.
- **The admissibility ledger pattern.** `claim_guards.py` / `claim_evidence.json` already model
  "visible failures, not just visible passes" — stage 7's output should be a ledger in the same
  spirit: every candidate that entered, where it died, why.

---

## 5. What's genuinely new and needs building, in build order

1. **Candidate-space enumeration + coarse pre-filter (stage 0–1).** Needs a decision on scope for the
   pilot: which axis set (current 256? the causal-attribution basis once the rotation test lands?),
   which target universe (curated pathways first, since that's cheapest and already has infra; genes
   individually is a much larger space and should come after the pilot validates the funnel shape),
   which cancer strata.
2. **The multiplicity-correction framework (stage 3), including the matched-random-control-at-capacity
   check from §2.1.** A predeclared search-space size, BH-FDR implementation reusing
   `honest_metrics`/`calibration` where possible, and a report format that states
   tested/passed-FDR/replicated counts the way `p2_floor_audit.py` reports fail/clear/unjudgeable —
   visible, not just the survivors.
3. **Novelty-filter feasibility scoping (§3, tier 2 and 3) — do this before promising it works.**
   Concretely: can PubMed E-utilities and Open Targets' API actually be reached from the environments
   this runs in, what are the rate limits, and what does a "no obvious prior hit" response actually
   mean in terms of false-negative rate (i.e., how often would a known association slip past the
   filter)? This is a scoping task, not a build task, and it should happen first because if it's not
   feasible, the whole "novel" claim needs a different, more honest framing (e.g., "not present in our
   own training targets" rather than "not published anywhere").
4. **A composition-attribution control (stage 5, per §1.1)** — reuse whatever cell-composition/purity
   baseline P1's own undischarged `composition_attribution` blocker calls for (ABSOLUTE purity is
   already wired for a related purpose in `d3_purity_result`), scored the same way any other confound
   control on this project is: must-fail if the candidate is composition-driven.
5. **The perturbation-consistency cross-check as a reusable pipeline step (stage 6)**, generalizing
   the bespoke analysis already done for E0/PBS into something that runs per-candidate rather than as
   one-off research code — scoped explicitly to cell-autonomous candidates (§1.1), not applied to
   composition/immune ones.
6. **The auto-drafted experiment proposal (stage 7).** Minimal viable version: a template that states
   the gene/pathway, the predicted direction and rough magnitude, a standard perturbation approach
   (CRISPR knockout/knockdown, matched to what the Perturb-seq resource already used, for
   consistency), a plausible cell line or model system for the cancer type in question, and the
   specific readout that would confirm or refute it. This is a design document per hypothesis, not
   something this project runs — the point is to hand a wet lab (or a collaborator) something
   immediately actionable rather than a bare correlation.

---

## 6. First concrete step — a small pilot, not mass scale on day one

Consistent with this project's own rule (`PROJECT_GUIDE.md` §2 rule 3: push results until they break,
don't accept the first favourable-looking thing) — but the same logic runs in reverse for a new
pipeline: **validate the funnel's shape at small scale before spending compute on mass scale.**

Recommended pilot, sized to be CPU/light-GPU feasible and not contend with training already in
flight:
- **Candidate space:** the curated Hallmark/KEGG pathway targets already in `frozen_rna_targets.npz`
  (known target space, no new data needed) crossed with the ~4–6 largest cancer strata, on the
  representation's existing 256 axes. Order of hundreds to low-thousands of cells — enough to
  exercise stages 0–4 for real, small enough to run today.
- **Explicitly skip stage 5's tier 2/3 (external novelty) in the pilot**, or run only tier 1
  (already-a-training-target) plus the composition-attribution control from §5.4 — the harder tiers
  need the feasibility scoping in §5.3 first, and running the pilot without it would produce
  hypotheses this plan cannot yet honestly call novel.
- **Deliverable:** a ledger of every candidate cell tested, how many cleared FDR (against both the
  permutation floor and the matched-random-control-at-capacity check from §2.1), how many replicated,
  how many survive the composition-attribution control, and — for whichever handful survive all of
  that — a first-draft stage-6/7 pass, to validate that the mechanistic cross-check and the
  experiment-proposal template actually produce something useful, not just that the statistics run.

An agent is running the §5.3 novelty-filter feasibility scoping and a first draft of the stage 0–4
pilot pipeline now (`NOTEBOOK_ENTRIES/` will carry the predeclaration and result, per standing
convention). Its scope predates §1.1/§2.1's findings and should be revisited against them when it
reports.

---

## 7. Honest risks

- **Tier-3 novelty filtering (literature search) cannot be made airtight**, and every hypothesis this
  pipeline surfaces must carry that caveat explicitly, not as a footnote once.
- **The multiple-testing burden at genuine mass scale (individual genes, not curated pathways) is
  large enough that FDR-controlled discovery may leave very few survivors** — this is not a flaw to
  route around, it's what an honest screen looks like. Report the attrition at every stage rather than
  only the final shortlist, the same discipline the rest of this project already uses.
- **The dominant morphological signal is very likely composition, not novel cell-autonomous biology
  (§1.1)** — expect most raw candidates to die at the composition-attribution control, and treat a
  funnel that doesn't lose most of its candidates there as a sign the control isn't working, not as a
  good yield.
- **The causal-grounding step (stage 6) only covers a fraction of what stage 2–4 can certify (§1.1)** —
  a statistically real, replicated, composition-cleared candidate outside the K562/RPE1 dictionary's
  domain is still reportable, but must be labelled as mechanistically ungrounded by this project's
  current resources, not silently dropped or silently promoted past what stage 6 actually showed.
- **A statistically real, replicated, mechanistically-consistent, apparently-novel association is
  still not a validated biological finding.** The furthest this pipeline can honestly go is a strong,
  specific, falsifiable hypothesis with a proposed experiment — the "genuine discovery" claim belongs
  to whoever runs that experiment, not to this project. Every output should say so.
