# Branch A — Literature Review and Chosen Direction

**MORPHEUS V2 · WSI→molecular-programme (Hallmark) prompting**
Master synthesis over 8 research lanes (A1_ledger) + adversarial novelty verdicts + real seed-42 artifacts.

---

## 1. Executive answer: is this a real, novel, usable contribution?

**PARTLY — yes as a *packaged mechanistic diagnosis*, no as any single headline claim.**

The honest verdict after adversarial adjudication is that **every individually-stated novelty claim that we tried to sell as a first-discovery was refuted**, but **the specific composite — a dual-head rank-collapse fingerprint tied to a control-adjusted, within-cancer Hallmark-prompting decomposition — survives** because no cited source assembles those pieces on this modality. Three of six adjudicated claims lost (2 refuted votes each); three survived. The contribution is real but must be pitched narrowly and defensively; leading with "we discovered WSI→molecular is confounded" gets sunk on sight by Buyer Beware (Nat. Biomed. Eng. 2026), DECAT (arXiv 2605.31504), and Howard 2021.

### Refuted novelty claims (do NOT claim these as ours)

| Claim | Lane | Why it fell | Cite instead |
|---|---|---|---|
| "Quantified confound decomposition (~50% cross-cancer) is novel" | a_prior_benchmarks | 2 refuted. Fu et al. 2020 and US Patent 12,299,884 already report the **same ~halving** (pan-cancer Spearman 62.8% → within-cancer 31.8%) and attribute it to cancer-type identity. | Fu 2020; Patent 12,299,884; MDPI Genes 2026; Buyer Beware 2026 |
| "Random-gene-set null control for WSI→Hallmark is novel" | c_random_control | 2 refuted. **Schmauch HE2RNA (2020) already ran 10,000 size-matched random gene lists as the specificity null in the exact WSI→Hallmark modality** and found real signatures beat random in only 25–50% of cancers. Venet 2011 / Domany 2018 own the general control. | Schmauch HE2RNA 2020; Venet 2011; Domany 2018; STAR Protocols 2022 |
| "Multi-objective head collapse onto a low-rank target is our phenomenon" | d_head_collapse | 2 refuted. NRC1 (NeurIPS 2024), Deep NRC (CPAL 2026, rank→intrinsic dim of low-rank target), Neural Collapse in Multi-Task Learning (OpenReview M4t2JUMlfI), and gradient-collapse-into-shared-subspace work cover both the single- and multi-objective rank-collapse framing. | Andriopoulos 2024; Rangamani & Unal 2026; Jing 2022; VICReg 2022 |

### Surviving novelty claims (these are OURS to keep)

| Claim | Lane | Status |
|---|---|---|
| **Dual-head mechanistic fingerprint**: biology rank ~5–6/256 collapsed while identity stays ~84, anchor residual ~0, as the *internal signature* of the cohort confound in a deployed WSI→molecular model. | b_batch_confound | **SURVIVES** (1 refuted, held). Matched-pair, within-model, per-objective collapse contrast is unreported. |
| **Packaged control-adjusted headline metric for WSI→Hallmark prompting**: within-cancer stratification × random-gene-set null → genuine specificity ~+0.07 that is *method-invariant* across all variants incl. baseline (SigLIP +0.005 the only residual). | e_benchmark_packaging | **SURVIVES** (0 refuted). The *combination* into a single benchmark metric is absent from all sources. |
| **Sharp empirical instantiation of the retrieval-vs-regression tradeoff** in WSI→molecular: matched dual-head minimal-sufficient collapse + the concrete diagnosis that a *per-dimension* VICReg variance floor without a covariance/decorrelation term cannot prevent alignment-induced rank collapse. | f_theory | **SURVIVES** (1 refuted, held). General theory is Wang CVPR'22 + FactorCL (cite, don't claim); the WSI instantiation + VICReg-mechanism diagnosis are ours. |

**One-line executive answer:** *Yes, there is a defensible, usable contribution — but it is the coupled "mechanism + benchmark" pairing (rank-collapse fingerprint ⟷ control-adjusted method-invariant specificity floor), not the confound, not the random control, and not collapse-as-phenomenon, each of which is prior art.*

---

## 2. Structured literature review by the 8 lanes

### Lane a — Prior WSI→molecular / pathway benchmarks
The canonical framework is **HEST-1k / HEST-Benchmark** (Jaume et al., NeurIPS 2024), which *structurally sidesteps* the cross-cancer confound by running one cancer type per task (patient-stratified, Pearson across patients) but never quantifies the pooled-vs-stratified gap. **PMC11814321** (Nat. Commun. 2025, 11 methods) documents low PCC (~0.28) and generalization failure, and explicitly notes it runs **no** mean-expression, random-gene, or permutation nulls — an acknowledged gap. **Arslan/Kather (Commun. Med. 2024)** report per-cancer AUC spread without isolating cancer-type-as-proxy. **MDPI Genes 2026** gives the closest qualitative statement ("high-performing models were capturing tumor type… as a proxy"). Critically, **Fu et al. 2020 (Nat. Cancer)** and **US Patent 12,299,884** already report the numeric ~50% halving (pan-cancer 62.8% → within-cancer 31.8%) — this is what refutes our decomposition-as-novel claim.

### Lane b — Batch / cohort confound in comp-path
**Howard et al. 2021 (Nat. Commun.)** is the canonical site-signature critique: preserved-site CV drops AUROC by mean 0.069 (max 0.291), 91% of features degrade, 36% vanish; site is predictable >0.85 AUROC after stain norm. **Dehkharghanian et al. 2023** recovers TCGA acquisition site at >86% from a cancer-type-trained network. **HESCAPE (arXiv 2508.01490)** finds a plain image encoder (Gigapath PCC 0.338) *matches or beats* aligned counterparts and names batch effects as dominant — its 0.338 baseline independently corroborates our global 0.348. **CHRep (arXiv 2604.21573)** names "regression-driven over-smoothing that suppresses biologically meaningful variation." **DECAT (arXiv 2605.31504)** is the most dangerous neighbour: CLIP embeddings falsely claim shared biology in 92% of confounded TCGA cases with within-cancer AUROC collapse (−0.24 TMB). None report our dual-head effective-rank asymmetry.

### Lane c — Random / null-model controls
Random / size-matched gene-set nulls are **codified standard practice**: AddModuleScore (Seurat) subtracts a size/expression-matched control set by construction; STAR Protocols 2022, NARGAB 2024, and eLife 2022 use random sets to measure specificity/false-positive rate. **Venet 2011** is canonical ("most random signatures are significantly prognostic"); **Domany 2018** extends it to 24/34 TCGA cancers and attributes the bias to cohort sub-classification (the transcriptomic analogue of our cross-cancer confound). **Ahlmann-Eltze 2025 (Nat. Methods)** — a mean baseline beats GEARS/scGPT, pretrained ≈ random-init — is the direct analogue of our "every method ~+0.07 incl. baseline." **Schmauch HE2RNA (2020)** already applied 10,000 random gene lists in the WSI→Hallmark modality — the finding that refuted our control-application-as-novel claim.

### Lane d — Multi-objective / dimensional head collapse
The single-objective core is a **named, citable phenomenon**: **NRC1** (Andriopoulos et al., NeurIPS 2024) — a d-dim representation regressed onto an n-dim target collapses to an n-dim subspace; **Deep NRC** (Rangamani & Unal, CPAL 2026) refines this to collapse toward the *intrinsic* dimension of low-rank targets (explaining rank 5–6 << nominal 50); **dimensional collapse** (Jing et al., ICLR 2022); **low-rank simplicity bias** (Huh et al.). Multi-task versions exist too — **Neural Collapse in Multi-Task Learning** (OpenReview M4t2JUMlfI) and gradient-collapse-into-task-agnostic-subspace work (arXiv 2509.16959) — which is why "multi-task interference is framed as accuracy not rank" was refuted. **VICReg** (Bardes et al., ICLR 2022) is explicit that the *per-dimension variance term does not prevent collapse* — the covariance/decorrelation term does. This is exactly why our 0.01 per-dim floor failed.

### Lane e — Benchmark packaging for adoption
The adopted-benchmark stack is a small repeatable set of layers: programmatic loader with fixed named splits (TDC `get_split(...)`), a standardized evaluator, a public leaderboard, and a machine-readable dataset card (**Croissant**, MLCommons Mar 2024, NeurIPS D&B recommended). Held-out test servers via containerized submission (MICCAI/Grand Challenge) or the HF private-evaluator four-object pattern are the leakage guard. **MorphoHELM (arXiv 2605.15383)** is the closest confound-aware template — hierarchical batch-stringency tiers with per-tier reporting — structurally identical to our global → within-cancer → random-control ladder, but on Cell Painting, not WSI→Hallmark. Current WSI→gene-expression benchmarks (biorxiv 2026.03.02.709012; HESCAPE) report raw per-gene correlation and **do not package** within-cancer stratification or a random-gene control. **This packaging gap is what our lane-e claim fills, and it drew zero refuted votes.**

### Lane f — Why one shared contrastive latent can't be optimal for retrieval + regression
The general proposition is prior art and should be cited, not claimed: **Wang et al. CVPR 2022** Theorem 1 (alignment drives toward the *minimal sufficient* representation, provably discarding non-shared task info I(v1,T|v2)); **FactorCL** (NeurIPS 2023, shared-vs-unique factorization); **Wang & Isola 2020** (alignment-uniformity tradeoff); **Jing 2022** + **LDReg** (ICLR 2024, effective rank bounds downstream capacity). Mapping to us: WSI = v1, Hallmark prompt = v2; the biology head aligned to the shared manifold and discarded the WSI-resident molecular variance not predictable from the coarse prompt — the +0.07 ceiling. The refuting note: this is a corollary of Theorem 1 + FactorCL. What survives is the **matched dual-head measurement** and the **VICReg per-dim-vs-covariance diagnosis**.

### Lane g — Venue / novelty scan
Closest three papers: **Buyer Beware / Nat. Biomed. Eng. 2026** (confounding critique legitimized at a top venue, but confounder axis is biomarker co-occurrence, not cancer-type Pearson decomposition, and no random-gene null); **HESCAPE** (alignment degrades expression prediction, batch named — but no Pearson decomposition, no rank analysis); **SupCon class-collapse theory (arXiv 2503.08203)** (predicts our mechanism but no empirical WSI instantiation). Reviewer verdict: novel as a *package* if positioned as a quantitative, control-adjusted, mechanistically-explained refinement — NOT as first-to-notice-confounding. Best-fit venues: **NeurIPS D&B** or **MICCAI** (benchmark + negative-control framing); **Nature Communications / NBME** if the mechanism is paired in; **ICLR** if rank-collapse is the headline.

### Lane h — Steelman the null
The strongest case against us: (i) the confound is named prior art (Buyer Beware, Howard) *with the fix already known* (stratified/site-preserved splits); (ii) **DECAT** already shows CLIP manufactures false shared biology + within-cancer collapse on TCGA — the most dangerous overlap; (iii) rank collapse is a textbook, fixable bug (mis-tuned VICReg floor) not a scientific finding; (iv) ~0.07 residual can be spun to declare the whole task uninteresting; (v) the SigLIP +0.005 edge is inside the noise band. **Residual novelty the null cannot erase:** the exact decomposed magnitudes on WSI→Hallmark *prompting* specifically, and the *mechanistic coupling* (alignment + weak per-dim floor → biology rank 5–6, identity healthy 84, anchor residual ~0) — which no single source ties together.

---

## 3. Gap map: what is genuinely unclaimed

Everything below is what remains after subtracting all cited prior art.

1. **[STRONGEST] The matched dual-head, per-objective rank-collapse contrast inside one deployed model.** Biology head rank ~5–6/256 collapsed vs sibling identity head healthy ~84, anchor residual ~0 (identity == frozen MLP-CLIP teacher). The collapse literature studies single objectives in isolation; nobody reports "one head collapses, its sibling doesn't, attributable to the auxiliary objective family (neighbour-KL + supcon)" as a within-model diagnostic. **Unclaimed.**

2. **[STRONG] The coupling of that fingerprint to benchmark method-invariance.** The claim that a rank-5–6 head still "works" on the benchmark *only because* ~46–49% of the signal is cross-cancer structure and random-gene controls already reach ~0.30–0.32, leaving ~+0.07 identical for every method. Prior art has the confound OR the collapse, never the causal bridge between them. **Unclaimed.**

3. **[MODERATE] The packaged control-adjusted headline metric for WSI→Hallmark *prompting*.** Within-cancer × random-gene-set null as the *primary reported number* (not raw global Pearson), on molecular-programme prompting. MorphoHELM has the template (other modality); HE2RNA has the random null (not packaged as a benchmark metric, not within-cancer-stratified, not prompting). The specific assembly is unclaimed. **Unclaimed.**

4. **[MODERATE] The VICReg per-dimension-vs-covariance failure diagnosis in a WSI→molecular alignment head.** VICReg theory predicts it; nobody has demonstrated the specific negative result — that a 0.01 per-dim variance floor fails to prevent alignment-induced rank collapse in this deployed setting. **Unclaimed as an empirical instantiation.**

**Explicitly NOT unclaimed (do not fight for these):** the existence of the cross-cancer confound; the ~50% halving magnitude; random-gene-set nulls as a control; NRC/dimensional collapse as a phenomenon; the alignment→minimal-sufficient theory; SigLIP-beats-CLIP.

---

## 4. Chosen formalization direction

**A METHOD/DIAGNOSTIC paper, anchored on the mechanism, using the benchmark decomposition as supporting evidence — NOT a pure benchmark paper and NOT a theory paper.**

### Justification (decisive)

- **A pure benchmark paper loses.** Lanes a, c, e's benchmark-novelty claims are the ones that got refuted or are weakest against DECAT/Buyer Beware/HE2RNA. A benchmark whose headline is "the task is confounded, here is the number" walks straight into the steelman (lane h) and the two refuted decomposition/control claims. The packaging (lane e) survived with 0 votes, but a benchmark alone is thin and adoption-gated (needs held-out server, Croissant, leaderboard — heavy lift for a modest ~+0.07 headline).

- **A pure theory paper loses.** Lane f is explicit: the proposition is a *corollary* of Wang CVPR'22 Theorem 1 + FactorCL. We cannot claim the theorem. A theory framing invites "already proven."

- **The method/diagnostic framing wins** because it is built on the two claims that *survived adjudication with the strongest residual novelty*: the dual-head rank-collapse fingerprint (lane b) and the VICReg-mechanism diagnosis (lane f), *coupled to* the method-invariance result (lane e). The mechanistic coupling is the one thing lane h (steelman) conceded it "cannot erase." It reframes the ~+0.07 result from "depressing negative benchmark" into "here is the internal fingerprint that predicts when your molecular-alignment head is riding cohort structure, and here is the specific regularizer term (covariance, not per-dim variance) that is missing." That is a *usable, prescriptive, mechanism-first* contribution the confound-critique papers (which are diagnosis-by-stratification only) do not offer.

**Direction, stated once:** *A diagnostic-and-mechanism paper: "effective-rank collapse of the biology head is the internal fingerprint of cohort-confounded WSI→molecular alignment," validated by a matched dual-head contrast, explained via minimal-sufficient/VICReg theory, and quantified against a random-gene-set-adjusted within-cancer specificity floor that is method-invariant.*

---

## 5. Concrete paper outline

**Working title:** *Rank Collapse as a Fingerprint of Cohort Confounding in WSI→Molecular-Programme Alignment*
(alt: *When the Biology Head Collapses: A Mechanistic Diagnostic for Confounded Histology-to-Hallmark Models*)

### Core claims (3)

**C1 (mechanism / headline).** In a shared multi-objective alignment model, the biology head that is aligned/regressed to a low-rank Hallmark manifold via neighbour-KL + supcon undergoes effective-rank collapse (rank ~5–6 of 256), while a sibling identity head with a rich frozen-teacher target stays healthy (rank ~84, anchor residual ~0). This matched within-model contrast is the internal fingerprint of the failure. *(Cite NRC1, Deep NRC, dimensional collapse for the phenomenon; the matched dual-head contrast is ours.)*

**C2 (coupling / benchmark consequence).** A biology head collapsed to rank ~5–6 still scores respectably on WSI→Hallmark prompting only because the benchmark is confounded: ~46–49% of global Pearson is cross-cancer cohort structure (0.348 → 0.188 within-cancer), a random-gene-set null already reaches ~0.30–0.32 globally, and the genuine within-cancer, random-control-adjusted specificity is only ~+0.07 — *identical across every method including the baseline*. *(Cite Buyer Beware, DECAT, Fu 2020, HE2RNA for the confound/control; the coupling of rank to method-invariance is ours.)*

**C3 (prescriptive diagnosis).** The collapse is theoretically expected (Wang CVPR'22 minimal-sufficient) and the mitigation failure is diagnosable: a *per-dimension* VICReg variance floor (0.01) provably cannot prevent it because rank is governed by the *covariance/decorrelation* term. *(Cite VICReg; the WSI-head empirical instantiation is ours.)*

*(Optional C4, weak — keep only as a minor result:* SigLIP is the only variant beating MLP-CLIP under the corrected metric (+0.005 within-cancer, wins 62% of targets); consistent with sigmoid-loss finer discrimination but inside the noise band.)*

### Experiments / tables and which of OUR results support each

| Table / Experiment | Content | Supporting OUR result |
|---|---|---|
| **T1 — Dual-head rank spectrum** | Effective rank of biology vs identity head across seeds; anchor-residual sanity check | **Biology rank ~5–6/256; identity ~84; anchor residual ~0 (identity == frozen MLP-CLIP teacher)** |
| **T2 — Confound decomposition ladder** | Global Pearson → within-cancer → random-gene-set control, per method | **Global 0.348 → within-cancer 0.188 (MLP-CLIP); ~46–49% cross-cancer; random control ~0.30–0.32** |
| **T3 — Method-invariance under adjustment** | Random-control-adjusted within-cancer specificity for all variants + baseline | **~+0.07 Pearson for EVERY method including baseline; SigLIP +0.005 residual, 62% target win** |
| **T4 — VICReg ablation** | Rank & specificity with per-dim variance floor on/off vs a covariance/decorrelation term added | **Per-dim variance floor (0.01) did NOT prevent rank collapse** (mechanism); covariance term as the predicted fix |
| **F1 — Fingerprint schematic** | 256-D biology head → ~50-D Hallmark manifold → intrinsic-dim collapse to 5–6, identity head preserved | Mechanism narrative tying C1→C2 |

### Claim → evidence mapping (explicit)

- **C1** ← T1 (rank ~5–6 vs ~84; anchor residual ~0), F1.
- **C2** ← T2 (0.348→0.188, ~46–49%, random ~0.30–0.32) + T3 (~+0.07 method-invariant).
- **C3** ← T4 (per-dim floor fails) + Wang/VICReg citations.
- **C4 (minor)** ← T3 (SigLIP +0.005, 62%).

### Minimal packaging (borrow from lane e, keep light)

Release the within-cancer split + **mandatory random-gene-set control run** + a Croissant dataset card, with the leaderboard/README reporting the **control-adjusted within-cancer delta** as the primary number. Do NOT gate the paper on building a held-out server — ship the loader + evaluator so the metric is reproducible, and position the benchmark as *supporting infrastructure for the diagnostic*, not the headline.

---

## 6. Honest risks (from the steelman-the-null lane) and defenses

| Risk (lane h) | Severity | Defense |
|---|---|---|
| **DECAT (arXiv 2605.31504) already shows CLIP manufactures false shared biology + within-cancer collapse on TCGA.** Adversary: "you re-derived DECAT on Hallmark labels." | **HIGHEST** | Cite DECAT prominently as concurrent/prior confirmation of the *confound half*. Differentiate on the *mechanism half*: DECAT works on classification AUROC over TMB/Age/TP53 and does not report a **matched dual-head effective-rank contrast**, a random-gene-set-adjusted **Pearson** floor for **molecular-programme prompting**, or the **VICReg per-dim-vs-covariance** diagnosis. Our headline is the rank fingerprint + prescriptive fix, not the confound. |
| **"Confound is known (Buyer Beware, Howard) and the fix (stratified splits) is standard — you rediscovered why to stratify."** | **HIGH** | Concede fully; never frame the confound as our discovery. Lead with C1 (mechanism) and C3 (fix), using C2 as corroboration. Emphasize we go *beyond* stratification: the rank fingerprint predicts the confound *from the model's internals* without needing to know the confounder a priori. |
| **"Random-gene control is standard (Venet, HE2RNA), and HE2RNA already ran it on WSI→Hallmark."** | **HIGH** | Cite HE2RNA and Venet/Domany as the origin of the control. Our claim is the *packaged within-cancer × random-null method-invariance* result, not the control's invention. |
| **"Rank collapse is a textbook, fixable bug (mis-tuned 0.01 floor) — engineering, not science."** | **MEDIUM-HIGH** | Turn it into the contribution: yes it is textbook (cite NRC/VICReg), and *that is the point* — we show the textbook failure fires in a deployed WSI→molecular head and, via T4, that the naive per-dim fix is the wrong knob. Frame as a diagnostic + correct-regularizer result, not a discovery of collapse. |
| **"~0.07 residual means the whole task is uninteresting — unpublishable negative."** | **MEDIUM** | Reframe: the value is the *method* (fingerprint + control-adjusted metric) that lets the field detect this on *any* WSI→molecular model, plus the mechanistic explanation for *why* method choice stops mattering. Method-invariance is itself a strong, falsifiable scientific claim. |
| **"Multi-task rank-collapse framing is also prior art (OpenReview M4t2JUMlfI, gradient-collapse work)."** | **MEDIUM** | Concede and cite. Keep the novelty strictly on the *matched within-model per-head contrast with anchor-residual control* in a deployed pathology model — not on multi-task collapse in the abstract. |
| **"SigLIP +0.005 is inside the noise band."** | **LOW** | Demote C4 to a minor result; do not hinge the paper on it. Report with the benchmark's ~0.2 variability caveat. |

---

## Recommended direction (≈200-word summary)

Write a **method/diagnostic paper**, not a benchmark or theory paper. The headline is a **mechanism**: in a multi-objective WSI→molecular alignment model, the biology head aligned to a ~50-D Hallmark manifold via neighbour-KL + supcon collapses to effective rank ~5–6 of 256, while a sibling identity head anchored to a frozen MLP-CLIP teacher stays healthy at ~84 (anchor residual ~0). This **matched, within-model, per-objective rank-collapse contrast** is the paper's novel core — the one thing the adversarial steelman conceded it cannot erase — and it serves as an *internal fingerprint* of cohort-confounded molecular alignment. Couple it to the benchmark consequence: the collapsed head still "works" only because ~46–49% of global Pearson is cross-cancer structure (0.348→0.188), a random-gene-set null already reaches ~0.30–0.32, and genuine within-cancer specificity is a method-invariant ~+0.07. Explain both via minimal-sufficient-representation theory (Wang CVPR'22, cite — don't claim) and diagnose the mitigation failure: a per-dimension VICReg floor cannot prevent collapse; the covariance/decorrelation term can. Concede the confound, the random control, and collapse-as-phenomenon to prior art (Buyer Beware, DECAT, HE2RNA, NRC1, VICReg). Ship a light within-cancer + random-null evaluator as supporting infrastructure. Target NeurIPS D&B or MICCAI; Nature Communications if the mechanism is foregrounded.
