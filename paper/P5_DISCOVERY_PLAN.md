# P5 — a discovery engine, not a query system: generating novel, testable biological hypotheses at scale

**Compiled:** 2026-08-05 · **Branch:** `research/rebase-vision` · **Status:** planning document, no
measurement yet. Written in response to a direct instruction: build a plan for genuine in-silico
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

## 1. The actual hard problem: "novel" is the part nothing on this project does yet

Everything else below — certified association, replication, mechanistic grounding — already has
working infrastructure (§3). **Filtering out rediscovery of known biology does not**, and it's the
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
any claim is made that the pipeline can tell novel from known — see §5.

---

## 2. The discovery funnel — mirrors a wet-lab screening funnel on purpose

Mass scale in, a short, defensible list out. Each stage is a real filter, and — matching this
project's ledger discipline — **every candidate that enters the funnel is tracked through it**,
including the ones that die at each stage, not just the survivors.

| stage | question | cost | reuses |
|---|---|---|---|
| **0 — enumerate** | what's the full candidate space? | cheap, combinatorial | — |
| **1 — coarse pre-filter** | cheap correlation scan, uncorrected, purely to cut compute | cheap | `calibra.spectral` |
| **2 — certify** | confound-adjusted, permutation null, injection-certified floor, **inductive** (out-of-sample, not transductive — see `PROJECT_GUIDE.md` §2 rule 3) | the existing instrument, at scale | `calibra.*` unchanged |
| **3 — correct for multiplicity** | BH-FDR across the *whole tested space*, predeclared before running | new — nothing on this project has swept an open space before; §4/§90-target work tested pre-chosen targets, not this | — |
| **4 — replicate** | does it hold in an independent cohort (ALCHEMIST; held-out TCGA cancer types)? | the existing infra | `alchemist_*`, `leave_sites_out.py` |
| **5 — novelty filter** | is this already known, at the three tiers in §1? | new, largely unbuilt | — |
| **6 — mechanistic ground** | does an independent functional-genomics resource predict the same direction? does a P3-certified causal name attach? | this project's specific asset — most projects doing this kind of screen don't have free perturbation data sitting on disk already | K562/RPE1 Perturb-seq loaders, `causal_attribution.CERTIFICATE` |
| **7 — score and propose** | rank by effect size over floor, replication strength, mechanistic consistency, novelty confidence, and **experimental feasibility** (standard cell line exists? gene is perturbable with standard tools?); auto-draft a minimal falsifiable experiment with a stated quantitative prediction | new | — |

**Stage 3 is the statistical crux and deserves its own paragraph.** Every certified result on this
project so far has tested a small, pre-chosen, frozen target list (~90 targets) — the multiple-testing
burden is bounded and known. An open sweep across genes/pathways/cancer-strata is a different
regime — potentially thousands to low-millions of candidate cells before the coarse pre-filter. The
search space size must be **fixed and stated before stage 2 runs**, exactly like every other
predeclaration on this project, so the FDR correction is honest rather than computed against
whatever survived. Two-stage design (FDR-controlled discovery + independent-cohort replication,
mirroring how genomics already does this) is the right shape — a candidate needs to clear *both*,
not either.

---

## 3. What already exists and plugs in directly — this is not starting from zero

- **The instrument.** `v2/calibra/` — confound adjustment (transductive and, as of today, inductive),
  the confound certificate, injection-certified detection/transmission floors, permutation nulls. Zero
  new statistics needed for stages 1, 2, or 4.
- **A second cohort for replication, already wired.** ALCHEMIST (1,106 paired NSCLC patients) — the
  channel already replicates there in aggregate (R=1.110). Per-candidate replication reuses the same
  machinery.
- **A free functional-genomics cross-check nobody else screening this way would have.** The K562/RPE1
  Perturb-seq resource (`e0_basis_transfer.py`, `pbs.py` loaders) already on disk. For any candidate
  gene covered by it, stage 6 is a real, in-silico, no-new-data-needed functional plausibility check:
  does perturbing this gene in vitro move things in the direction the hypothesis predicts?
- **A mechanistic-naming certificate.** P3's `causal_attribution.CERTIFICATE` (four conditions,
  currently 29/128 axes certified, an agent testing right now whether a better basis raises that).
  Where it applies, a candidate gets a causal story, not just a correlation.
- **The admissibility ledger pattern.** `claim_guards.py` / `claim_evidence.json` already model
  "visible failures, not just visible passes" — stage 7's output should be a ledger in the same
  spirit: every candidate that entered, where it died, why.

---

## 4. What's genuinely new and needs building, in build order

1. **Candidate-space enumeration + coarse pre-filter (stage 0–1).** Needs a decision on scope for the
   pilot: which axis set (current 256? the causal-attribution basis once §3's rotation test lands?),
   which target universe (curated pathways first, since that's cheapest and already has infra; genes
   individually is a much larger space and should come after the pilot validates the funnel shape),
   which cancer strata.
2. **The multiplicity-correction framework (stage 3).** A predeclared search-space size, BH-FDR
   implementation reusing `honest_metrics`/`calibration` where possible, and a report format that
   states tested/passed-FDR/replicated counts the way `p2_floor_audit.py` reports fail/clear/
   unjudgeable — visible, not just the survivors.
3. **Novelty-filter feasibility scoping (§1, tier 2 and 3) — do this before promising it works.**
   Concretely: can PubMed E-utilities and Open Targets' API actually be reached from the environments
   this runs in, what are the rate limits, and what does a "no obvious prior hit" response actually
   mean in terms of false-negative rate (i.e., how often would a known association slip past the
   filter)? This is a scoping task, not a build task, and it should happen first because if it's not
   feasible, the whole "novel" claim needs a different, more honest framing (e.g., "not present in our
   own training targets" rather than "not published anywhere").
4. **The perturbation-consistency cross-check as a reusable pipeline step (stage 6)**, generalizing
   the bespoke analysis already done for E0/PBS into something that runs per-candidate rather than as
   one-off research code.
5. **The auto-drafted experiment proposal (stage 7).** Minimal viable version: a template that states
   the gene/pathway, the predicted direction and rough magnitude, a standard perturbation approach
   (CRISPR knockout/knockdown, matched to what the Perturb-seq resource already used, for
   consistency), a plausible cell line or model system for the cancer type in question, and the
   specific readout that would confirm or refute it. This is a design document per hypothesis, not
   something this project runs — the point is to hand a wet lab (or a collaborator) something
   immediately actionable rather than a bare correlation.

---

## 5. First concrete step — a small pilot, not mass scale on day one

Consistent with this project's own rule (§2 of `PROJECT_GUIDE.md`: push results until they break,
don't accept the first favourable-looking thing) — but the same logic runs in reverse for a new
pipeline: **validate the funnel's shape at small scale before spending compute on mass scale.**

Recommended pilot, sized to be CPU/light-GPU feasible and not contend with training already in
flight:
- **Candidate space:** the curated Hallmark/KEGG pathway targets already in `frozen_rna_targets.npz`
  (known target space, no new data needed) crossed with the ~4–6 largest cancer strata, on the
  representation's existing 256 axes. Order of hundreds to low-thousands of cells — enough to
  exercise stages 0–4 for real, small enough to run today.
- **Explicitly skip stage 5 (novelty) in the pilot**, or run only tier 1 (already-a-training-target)
  — the harder tiers need the feasibility scoping in §4.3 first, and running the pilot without it
  would produce hypotheses this plan cannot yet honestly call novel.
- **Deliverable:** a ledger of every candidate cell tested, how many cleared FDR, how many
  replicated, and — for whichever handful survive both — a first-draft stage-6/7 pass, to validate
  that the mechanistic cross-check and the experiment-proposal template actually produce something
  useful, not just that the statistics run.

An agent is being sent on the §4.3 novelty-filter feasibility scoping and a first draft of the
stage 0–4 pilot pipeline now (`NOTEBOOK_ENTRIES/` will carry the predeclaration and result, per
standing convention).

---

## 6. Honest risks

- **Tier-3 novelty filtering (literature search) cannot be made airtight**, and every hypothesis this
  pipeline surfaces must carry that caveat explicitly, not as a footnote once.
- **The multiple-testing burden at genuine mass scale (individual genes, not curated pathways) is
  large enough that FDR-controlled discovery may leave very few survivors** — this is not a flaw to
  route around, it's what an honest screen looks like. Report the attrition at every stage rather than
  only the final shortlist, the same discipline the rest of this project already uses.
- **A statistically real, replicated, mechanistically-consistent, apparently-novel association is
  still not a validated biological finding.** The furthest this pipeline can honestly go is a strong,
  specific, falsifiable hypothesis with a proposed experiment — the "genuine discovery" claim belongs
  to whoever runs that experiment, not to this project. Every output should say so.
