# MORPHEUS Rebase — Ranked Research Paths

*Decisive, honest recommendation. Every A1–A5 *mechanism* is pre-empted (see LITERATURE_REVIEW); all 8
adversarially-vetted candidates SURVIVED but each has a reframe-or-die near-collision (see NOVELTY_LEDGER).
The paths below are ranked by **defensibility × whitespace × MORPHEUS-asset support**, not by ambition.
The winning strategy is deliberate: lead with the two contributions that are hardest to scoop and that
directly repair the project's diagnosed wound — **evaluation** (A3) and the **unproven identifiability→reliability
link** (A2) — then let the promptable-interface (A1) and the system paper (A2×A4×A5) build on that spine.*

**MORPHEUS assets referenced throughout** (from REBASE_THESIS §2): `v2/TumorStateV2` (WSI+RNA aligned; typed
identity/biology/context/patient heads; `DenseAdapter`+availability-mask fusion; leakage-safe targets;
safety-gated closed-RAG cards) · PRISM `models/tqi.py` (unwired `ScopeDetector` + head-registry scaffold) ·
`QueryBlock` typed-slot attention · `text_prototypes` embedding infra · CPTAC proteomics/phospho inventory ·
the session's **rank-decoupling / covariance-term finding** (biology head collapses to low rank; a covariance
term reverses it +53 / ~2.1× over 3 seeds; confounded benchmark score is *unchanged* → the eval is structurally
blind to representation quality).

---

## Tier 1 — Recommended lead papers (build these first)

### PATH 1 (rank #1) — BioELK-Bench: a confound-aware emergent-biology elicitation benchmark with a validity certificate
*Ledger C7 (A3), with C8 (A3×A5) as its causal companion.*

**Claim.** Ship the evaluation instrument the field lacks: a curated probe battery measuring how much
pathway/mechanism knowledge a frozen multimodal tumor trunk encodes and how readily NL prompting elicits it,
scored to survive four killer critiques — (i) MDL/MI + thresholded accuracy under smooth AND hard metrics;
(ii) Hewitt–Liang selectivity vs a random-programme control; (iii) a Farquhar distractor/site sanity check
(the recovered direction must FAIL to predict submitting-site/scanner while succeeding on biology); (iv) a
decisive delta vs the GenePT text-prior and Kedzierska PCA baselines; (v) a password-locked control separating
"trunk lacks knowledge" from "prompt failed." Emergence plotted loss-indexed, at per-pathway quantum grain,
across multiple elicitation formats. **Companion (C8):** an amnesic-counterfactual *causal-use validity gate* —
erase a pathway direction (INLP/LEACE) and confirm the survival/drug-response prediction changes; use it as the
accept/reject gate on the closed-RAG hypothesis card, defusing the ELK human-simulator hazard.

**Novelty vs prior art (from ledger).** Every ingredient is borrowed (Hewitt–Liang, Voita–Titov MDL, Farquhar,
Greenblatt password-lock, Elazar amnesic, Med-PaLM rubric) but the *assembled* confound-certified,
text-prior-beating elicitation battery **for a multimodal cancer FM does not exist**. Closest: "Probing, Fusion,
and Trustworthiness" (2606.17115) probes a multimodal cancer FM but has no certificate/selectivity/MDL/NL-elicitation/
text-prior/password-lock/emergence-curves; VCBench is single-cell, not NL-elicitation. **The validity certificate
(probe reads biology, not the most-prominent site/distractor) is the unfilled contribution, and it ties A3
evaluation to A2 identifiability.**

**Required experiments.** (1) Assemble the probe battery: per-pathway readouts from curated gene-set/programme
labels, each paired with a random-programme control task. (2) Implement MDL + thresholded scoring; report both.
(3) Farquhar site/scanner sanity: train site-predictors on each recovered direction; require failure on site,
success on biology (report a de Jong Robustness-Index-style scalar). (4) Beat GenePT text-prior + HVG+PCA
decisively. (5) Password-lock control: fine-tune a locked capability, verify the battery distinguishes
un-elicited from absent. (6) Loss-indexed emergence curves across ≥3 elicitation formats (anti-sandbagging).
(7) C8 gate: LEACE/INLP erasure → measure Δ on survival/drug-response; steer to confirm sufficiency; wire as the
hypothesis-card accept/reject rule.

**MORPHEUS assets that support it.** The **rank-decoupling finding** is the motivating result — this benchmark is
exactly the instrument that would have *seen* the covariance-term rank recovery the confounded C-index missed.
`text_prototypes` supplies the NL elicitation formats; the **safety-gated closed-RAG cards** are the artifact the
C8 causal gate certifies; leakage-safe targets + the typed biology head supply probe labels.

**Build-on / new-design split.** *Build on:* CCS/difference-in-means/SAPLMA readouts, Hewitt–Liang/Voita–Titov/
Farquhar/Greenblatt/Elazar/LEACE methods, Med-PaLM multi-axis rubric, GenePT/Kedzierska nulls. *New design:* the
confound-**certificate** for a multimodal-cancer trunk; per-pathway-quantum loss-indexed emergence curves for
oncology; the amnesic causal-use gate on a generated hypothesis card.

**Why #1.** Highest defensibility (an evaluation contribution is hardest to scoop), clearest whitespace ("genuinely
unoccupied"), lowest dependence on unproven modeling claims, and it *repairs the diagnosed wound* — the harness
blind to representation quality. It is also the prerequisite instrument that makes Paths 2–4 falsifiable.

---

### PATH 2 (rank #2) — Identifiability-as-reliability: proving the missing empirical link
*Ledger C4 (A2). The load-bearing claim of the entire rebase.*

**Claim.** Demonstrate that *identified* pathway slots (iVAE-style conditional prior + sparse-mechanism-shift on
perturbation data) yield measurably MORE reliable and transferable NL prompting than (i) a-priori gene-set pathway
tokens (SurvPath) and (ii) unsupervised/linear factors (MOFA+/PCA), on **leave-site-out and leave-cancer-out**
splits under a confound-aware protocol. Ship the reliability instrument: reduction in Dawood-style biomarker
co-dependence leakage, de Jong Robustness-Index gain, Roeder cross-run linear-identifiability stability.

**Novelty vs prior art (from ledger).** Every prior (SurvPath, SENA, biolord, CPA, COMPASS) claims *interpretability*
or *accuracy*; **none demonstrates that identifiability itself buys prompt reliability/transfer.** Closest: SurvPath
(pathway tokens, but a-priori/fixed, claims interpretability not identifiability-driven reliability). This is the
**least pre-empted claim in the corpus** — the "identifiability→reliability" arrow is drawn by no one.

**Required experiments.** (1) Three matched slot representations on the same frozen WSI+RNA trunk: (a) identified
(conditional-prior + sparse-mechanism-shift), (b) SurvPath a-priori gene-set tokens, (c) MOFA+/PCA factors.
(2) Prompt each identically; measure reliability (prompt-paraphrase consistency, cross-run stability via Roeder
linear-identifiability) and transfer (leave-site-out, leave-cancer-out C-index/AUC delta). (3) Co-dependence-leakage:
Dawood-style test that a single-biomarker prompt returns the isolated programme, not the correlated bundle — does
identification reduce leakage? (4) Robustness-Index: does the identified latent raise RI>1 (biology-neighbours beat
same-center neighbours)? (5) Report perturbation coverage against the block-identifiability ceiling (JMLR-2024;
Squires 1/node). Claim only **block/programme-group** identifiability.

**MORPHEUS assets that support it.** The biology head + programme regression are the substrate to identify;
`DenseAdapter`+availability-mask fusion supplies the conditional-prior auxiliary variables (site/tissue/patient
metadata as iVAE `u`); the cancer-held-out guardrail is already the required split; the rank-decoupling finding is
the pilot evidence that representation quality varies invisibly to current scoring — this path *measures* it.

**Build-on / new-design split.** *Build on:* iVAE conditional-prior identification, sVAE+/SAMS sparse-mechanism-shift,
SENA pathway projection, Dawood/de Jong/Roeder instruments, SurvPath/MOFA+ baselines. *New design:* the **causal test
that identifiability (not interpretability) buys reliability/transfer**, and the reliability instrument itself.

**Why #2.** It is MORPHEUS's *defensible core* and directly answers the reviewer's "isn't this just SurvPath/SENA?".
Ranked below Path 1 only because it depends on Path 1's benchmark to score reliability credibly, and because its
modeling machinery (identification) is more scoopable than a pure evaluation contribution.

---

## Tier 2 — Strong follow-on papers (build on Tier 1)

### PATH 3 (rank #3) — Auto-detecting & abstaining on under-specified biological queries
*Ledger C1 (A1), with C2 (A1×A2) as the routing-mechanism companion.*

**Claim.** A router that infers the latent biological task from an *ambiguous* NL query over the frozen
WSI+molecular tumor state, returns *multiple valid programme interpretations* when under-specified (SAM-style),
and *abstains* when no identified pathway slot addresses it — formalized as a posterior over identified programmes
(Xie implicit-Bayes), evaluated on a **held-out-task-family** split. Companion C2: internalize the routing inside
the frozen trunk (Skill-MoE-style NL-skill routing to *named identified pathway slots*, task materialized as a
Hendel/Todd function vector), short-horizon (one prompt → one addressed programme).

**Novelty vs prior art (from ledger).** Every in-domain promptable system (Med-PaLM M, PathChat, ChatNT, LangCell)
is *told* the task and none abstains. Closest: Nguyen "Answer/Clarify/Abstain" med-VLM triage — but its latent is
*input-quality*, not *biological task*; radiology, not WSI+molecular; no multi-interpretation; not held-out-task-family.
**Differentiator: the detected latent is the biological task/programme, with multi-interpretation + abstention over
*identified* slots, evaluated across held-out task families.**

**Required experiments.** (1) Build/borrow a declarative biological task census (Ledger C3: hundreds of tumor-state
asks, many paraphrases, held-out families). (2) Router head over the identified slots (Path 2) emitting a posterior;
abstain when max-posterior < calibrated threshold; multi-output when the posterior is multi-modal. (3) Held-out-task-family
generalization (FLAN/Super-NaturalInstructions protocol) — the credibility test vs memorized heads. (4) Calibration
(Kadavath P(True)/P(IK)) and abstention accuracy-coverage curves. (5) C2: benchmark internalized routing consistency/
identifiability against an external agent toolbox (Biomni-style) — answers the "how is this not Biomni?" objection.

**MORPHEUS assets that support it.** PRISM `tqi.py` `ScopeDetector`+head-registry is the exact scaffold to wire;
`QueryBlock` typed-slot attention is the routing substrate; `text_prototypes` provides the NL front-end; the
identified slots come from Path 2.

**Build-on / new-design split.** *Build on:* Xie implicit-Bayes, SAM ambiguity-as-multi-output, InstructBLIP
query-conditioned readout, Hendel/Todd task/function vectors, FLAN held-out-family protocol, Skill-MoE NL routing.
*New design:* auto-detection + multi-interpretation + **abstention over biologically-identified pathway slots**;
internalized (vs agent-glue) routing to *named identified* slots.

**Why #3.** The sharpest A1 delta, and it is the natural *interface* over the Path-2 slots and the Path-1 benchmark.
Ranked below Tier 1 because it *depends* on identified slots existing (Path 2) and because its closest neighbour
(Nguyen triage) is a live, recent reframe risk requiring careful positioning.

---

### PATH 4 (rank #4) — Frozen-trunk, block-identified pathway slots: the system paper
*Ledger C5 (A2×A4×A5). The "build the thing" integration paper.*

**Claim.** A lightweight adapter on the FROZEN multimodal trunk (WSI+RNA + optional proteomics/phospho) that
identifies pathway-GROUP (block) latent slots via sparse-mechanism-shift (sVAE+/SAMS) + known pathway-membership
grouping (Morioka) + known pathway causal graph (CauCA) as the named Locatello bias, projects them onto named
pathways (SENA-style), and exposes them for free-form NL addressing — claiming only *partial/programme-group*
identifiability against the JMLR-2024 ceiling.

**Novelty vs prior art (from ledger).** Closest: Winter et al. 2606.29949 (frozen WSI+RNA + alignment adapter +
open-vocab NL pathway prompting + honest "graduated predictability spectrum") **anticipates the system combination
but not the identifiability guarantee**; SENA-discrepancy-VAE has the guarantee but is single-omic, non-frozen, no
NL adapter. **Differentiator: the block-identifiability *guarantee* on identified *latent* slots exposed through a
frozen-multimodal NL adapter — the union neither holds.**

**Required experiments.** (1) Adapter identifying block slots on the frozen trunk using CPTAC/perturbation signal +
pathway grouping + pathway graph. (2) SENA pathway projection for naming. (3) Expose slots to the Path-3 NL router.
(4) Prove block-identifiability empirically (cross-run stability, perturbation coverage vs Squires/Varıcı budgets).
(5) A4 sub-study: for proteomics/phospho, encode-vs-retrieve ablation (the un-voted MES direction — beat In-Context-RALM
zero-surgery + GenePT null, show gap survives distillation). (6) A5 sub-study: forward + inverse (PDGrapher-style)
counterfactuals over the slots with calibrated abstention, beating the linear mean-shift null (Virtual Cell Challenge
protocol).

**MORPHEUS assets that support it.** `DenseAdapter`+availability-mask fusion is the mosaic/frozen-plug-in substrate;
CPTAC proteomics/phospho inventory is the marquee encode-vs-RAG test modality; `QueryBlock` exposes `(batch, n_programme, D)`
slots the unwired `tqi.py` needs.

**Build-on / new-design split.** *Build on:* sVAE+/SAMS/SENA/CauCA/Morioka identification, BLIP-2 Q-Former plug-in,
MoPoE mosaic fusion, PDGrapher inverse queries, conformal/GPerturb UQ. *New design:* the identifiability *guarantee*
on the frozen-multimodal-NL substrate; the pathway-slot-keyed molecular-memory encode-vs-retrieve decision.

**Why #4.** The most ambitious and most useful as a flagship system, but ranked last in the recommended set because
its novelty is the *most pre-empted* (Winter + SENA squeeze it), it *depends* on Paths 1–3, and it folds in
un-adjudicated A4/A5 directions that should be voted first (NOVELTY_LEDGER §3).

---

## Tier 3 — Supporting / infrastructure (do, but not standalone flagship)

- **PATH 5 — Declarative biological task census + emergence-conditioned corpus (Ledger C3, A1×A3).** Ship as the
  *dataset/benchmark backbone* that Paths 1 and 3 consume (held-out task families; corpus curated for Chan et al.'s
  burstiness/long-tail preconditions; report auto-detection *as a function of* measured statistics). High value as
  infrastructure; weaker as a standalone headline. Closest: SC-Arena + Biology-Instructions.
- **PATH 6 — Context-conditioned compositional slots + biological-Homomorphism-Error objective + VSA capacity
  budget (Ledger C6, A2×A5).** A rigorous, measurement-forward A2 deepening (programme×cell-state slots, HE objective,
  measured binding ceiling via a ConceptMix-style combinatorial-load benchmark). Best as a follow-on once Path 2
  establishes identified slots; first *biological* HE instantiation is the hook.

---

## Explicitly retired (do NOT lead with these — see NOVELTY_LEDGER §4)

No adjudicated candidate was refuted, so nothing is *killed* — but the **bare axis framings** are dead as headlines
and must not anchor any paper:
- "An NL-promptable / task-inferring multimodal cancer model" (Med-PaLM M, PathChat, ChatNT).
- "Pathway-addressable multimodal tumor tokens" as *architecture* (SurvPath) — survives only as Path 2's
  *identified-buys-reliability* claim.
- "Perturbation/drug as a query on a learned latent" (scGen/CPA/GEARS/biolord/STATE/Geneformer) — survives only as
  Path 4's frozen/identified/measured-geometry substrate.
- "Encode-vs-retrieve mechanism / frozen-trunk-plus-datastore" (kNN-LM/RAG/In-Context-RALM) — survives only as the
  molecular-modality *decision rule* (un-voted; vote before committing).
- "Our model has emergent biological knowledge" (Schaeffer mirage; Kedzierska/Ahlmann-Eltze nulls) — survives only as
  Path 1's *measurement method with a validity certificate*.

## Sequencing recommendation

Build **Path 1 (benchmark)** and **Path 2 (identifiability→reliability)** in parallel — they are mutually reinforcing
and constitute the defensible spine. **Path 3 (auto-detect+abstain)** then layers the NL interface over Path 2's slots
and is scored by Path 1. **Path 4 (system paper)** integrates everything and should wait until the un-voted A4-MES and
A5-identification-buys-intervention directions clear their own adversarial vote. Paths 5–6 are infrastructure/deepening
that fall out of the first three.
