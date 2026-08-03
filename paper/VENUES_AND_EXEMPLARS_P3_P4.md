# Submission venues and structural exemplars — P3 and P4

**Compiled:** 2026-08-03 · **Branch:** `research/rebase-vision`

---

## 0. How to read this document (verification protocol)

Three fabricated citations have already contaminated this project. This document was
therefore built under a retrieval-only rule: **no citation here was written from memory.**

Every paper below was retrieved during compilation from one of:

- Crossref REST API (`https://api.crossref.org/works/<DOI>`) — authoritative title, author
  list, container title, volume, issue, pages, date, type;
- Europe PMC REST API (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`) — record
  metadata, open-access status, PMCID;
- Europe PMC full-text XML (`.../rest/<PMCID>/fullTextXML`) — **actual section heading
  order, figure/table/reference counts, body word counts, and all verbatim quotations**;
- arXiv API (`https://export.arxiv.org/api/query`) — preprint metadata and comments;
- direct fetch of publisher policy pages (URL given inline in each case).

Each claim is tagged:

- **[V]** = verified — read directly from a retrieved record, full text, or policy page.
- **[I]** = inferred — a reasonable reading not stated in the source. Treat as an opinion.
- **[U]** = unverified — could not be retrieved in this session. **Check before relying on it.**

Section 7 is a ledger of everything that stayed **[U]**, so the gaps are visible rather
than buried.

> A note on method: WebSearch was unavailable for this compilation, so discovery ran
> through structured bibliographic APIs instead of a search engine. This narrows
> serendipity but raises precision — every hit arrived as a metadata record with a DOI
> attached, so there was no step at which a plausible-sounding title could enter the
> document without a resolvable identifier behind it.

---

## 1. The decisive question for P3 — bad news first

**Asked:** is a well-controlled negative result on supervision-target choice publishable in
the *computational-pathology mainstream*, or does it need a venue with an explicit
negative-results policy?

**Answer: neither, and the framing of the question is the trap.**

Three findings, in descending order of how much they should change the plan.

### 1.1 The computational-pathology imaging mainstream has no pathway for this paper

I searched for published papers in the MICCAI / *Medical Image Analysis* orbit whose
**headline** is "X does not work better than Y". The search was deliberately generous
(Crossref filtered to `issn:1361-8415` for MedIA; arXiv full-text queries over
`eess.IV`; Europe PMC title queries for `"do not outperform"`, `"fail to outperform"`,
`"no better than"`). What came back:

- **MedIA does publish critical work**, but rarely and without a dedicated article type.
  The one clean example retrieved: Liu B, Dolz J, Galdran A, Kobbi R, Ben Ayed I. *Do we
  really need dice? The hidden region-size biases of segmentation losses.* **Medical Image
  Analysis** 91:103015 (2024). DOI [10.1016/j.media.2023.103015](https://doi.org/10.1016/j.media.2023.103015). **[V]**
- **MICCAI 2026's call for papers** (fetched from
  `https://conferences.miccai.org/2026/en/CALL-FOR-PAPERS.html`) is themed **"From
  Algorithms to Clinical Translation"** and solicits **"translational solutions to
  real-world clinical challenges"**. It lists **no separate tracks or article types**, and
  makes **no mention** of negative results, validation studies, benchmarking, or
  reproducibility as welcomed contributions. **[V]**
- The most on-point computational-pathology negative result I could find — Ganz J,
  Ammeling J, Rosbach E, Lausser L, Bertram CA, Breininger K, Aubreville M. *Is
  Self-supervision Enough? Benchmarking Foundation Models Against End-to-end Training for
  Mitotic Figure Classification.* In *Bildverarbeitung für die Medizin 2025*, Informatik
  aktuell, Springer Fachmedien Wiesbaden, pp. 63–68 (2025). DOI
  [10.1007/978-3-658-47422-5_15](https://doi.org/10.1007/978-3-658-47422-5_15) — landed in a
  **regional German workshop proceedings**, six pages long. **[V]** Its abstract states
  the finding plainly: *"We found that the end-to-end-trained baseline outperformed all
  FM-based classifiers, regardless of the amount of data provided… rendering both of the
  above assumptions incorrect."* (verbatim, arXiv:2412.06365 abstract) **[V]**

  The instructive part is what happened next. The **extended** version of that work, by
  largely the same group, appeared in *Machine Learning for Biomedical Imaging* vol. 3,
  MELBA–BVM 2025 Special Issue (2026), DOI
  [10.59275/j.melba.2026-a3eb](https://doi.org/10.59275/j.melba.2026-a3eb) — retitled
  *Benchmarking Foundation Models for Mitotic Figure Classification*, with a **positive**
  headline (LoRA-adapted foundation models perform well). **[V]** The negative framing got
  six pages in a workshop; the positive framing got the journal. Draw your own conclusion,
  but do not assume it was a coincidence.

**Implication for P3:** submitting P3 as a computational-pathology *methods* paper means
competing in a venue whose stated theme is clinical translation, in an 8-ish-page
conference format **[U — page limit not verified, see §7]**, against reviewers with no
article type to slot it into. The three-seed design, the bootstrap CIs, the negative
control, the matched-capacity PCA arm and the placebo deflation — the entire evidentiary
apparatus that makes P3 decisive — is exactly what does not fit in that format.

### 1.2 The good news: P3's genre is thriving, one field over

The "sophisticated method loses to a deliberately simple baseline, and the paper is
*about* that" genre has **strong, recent, high-prestige precedent in computational
biology and ML-for-biology** — not in medical imaging. All verified:

| Paper | Venue | Type | Verified via |
|---|---|---|---|
| Ahlmann-Eltze C, Huber W, Anders S. *Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.* Nat Methods **22**(8):1657–1661 (2025). DOI [10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6) | **Nature Methods** | **Brief Communication** (`<subject>` tag in JATS full text) | Crossref + PMC12328236 full text **[V]** |
| Kedzierska KZ, Crawford L, Amini AP, Lu AX. *Zero-shot evaluation reveals limitations of single-cell foundation models.* Genome Biol **26**(1):101 (2025). DOI [10.1186/s13059-025-03574-x](https://doi.org/10.1186/s13059-025-03574-x) | **Genome Biology** | Research | Crossref + PMC12007350 full text **[V]** |
| Venet D, Dumont JE, Detours V. *Most Random Gene Expression Signatures Are Significantly Associated with Breast Cancer Outcome.* PLoS Comput Biol **7**(10):e1002240 (2011). DOI [10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240) | **PLOS Comp Biol** | Research | Crossref + PMC3197658 full text **[V]** |
| Christodoulou E, Ma J, Collins GS, Steyerberg EW, Verbakel JY, Van Calster B. *A systematic review shows no performance benefit of machine learning over logistic regression for clinical prediction models.* J Clin Epidemiol **110**:12–22 (2019). DOI [10.1016/j.jclinepi.2019.02.004](https://doi.org/10.1016/j.jclinepi.2019.02.004) | **J Clin Epidemiol** | Systematic review | Crossref **[V]** |
| Raghu M, Zhang C, Kleinberg J, Bengio S. *Transfusion: Understanding Transfer Learning for Medical Imaging.* arXiv:1902.07208; arXiv comment field reads **"NeurIPS 2019"** | **NeurIPS** | Conference paper | arXiv API **[V]**; peer-reviewed status inferred from the comment field **[I]** |
| Boettcher S. *Inability of a graph neural network heuristic to outperform greedy algorithms in solving combinatorial optimization problems.* Nat Mach Intell **5**:24–25 (2022). DOI [10.1038/s42256-022-00587-0](https://doi.org/10.1038/s42256-022-00587-0) | **Nature Machine Intelligence** | 2 pages → almost certainly Matters Arising **[I]** | Crossref **[V]** |

So: **P3 does not need a "negative results" journal.** *Nature Methods* put this exact
genre in print in 2025 with the word "not" in the title. What P3 needs is to be
**submitted as a computational-biology benchmarking/analysis paper, not as a
computational-pathology methods paper.**

### 1.3 The uncomfortable structural point

P3's supervision target is a Perturb-seq-derived dictionary and its comparator is Hallmark
gene sets. Its evaluation is on held-out molecular targets. That is a **computational
biology** contribution wearing a histology encoder. The imaging community will read it as
"an ablation"; the comp-bio community will read it as "a benchmarking result about whether
interventional data buys you a better supervision signal" — which is a live, contested
question in that field right now, as §2.1 shows.

**Recommendation: reposition, don't downgrade.** Aim high in comp-bio, keep MELBA as the
comp-path-community fallback, and treat the explicit-negative-results venues (PLOS ONE,
GigaScience) purely as a floor.

### 1.4 One thing to prepare for

The Ahlmann-Eltze *Nature Methods* paper drew a **public, titled rebuttal**: Miller HE,
Mejia GM, Leblanc FJA, Wang B, Swain B, de Lima Camillo LP. *Deep Learning-Based Genetic
Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics.*
bioRxiv, 2025-10-21. DOI
[10.1101/2025.10.20.683304](https://doi.org/10.1101/2025.10.20.683304) **[V]**

That is the "you just did it wrong" objection, materialised, in public, against the best
paper in this genre — and its counter-argument is **about the choice of metric**. P3's
analogous exposure is the choice of evaluation targets and the matched-capacity
construction. Whatever the final framing, P3 should pre-register in its own text why the
40 held-out targets and the capacity matching are the right yardstick, and should say so
in the Results, not only the Methods. See §3.3 for how the exemplars do this.

---

## 2. P3 — ranked venues

### Rank 1 — *Genome Biology* (BMC/Springer Nature)

- **Scope fit:** high. Perturb-seq dictionaries, Hallmark gene sets, held-out molecular
  targets, random-gene-set negative controls — this is squarely the readership. **[I]**
- **Negative results:** **strong precedent, not policy.** Kedzierska et al. 2025 (Genome
  Biol 26:101) is a full research article whose title is *"Zero-shot evaluation reveals
  **limitations** of single-cell foundation models"* and whose abstract states the models
  *"could be outperformed by simpler methods"* (verbatim, PMC12007350 abstract). **[V]**
  No explicit negative-results *policy* was retrieved. **[U]**
- **Length:** the exemplar runs ~4,900 body words, 2 main figures, 45 references, with
  supplementary carrying the rest. **[V]** Formal limits not retrieved. **[U]**
- **Open access / APC:** the exemplar is open access. **[V]** APC amount not retrieved. **[U]**
- **Time to first decision:** not retrieved. **[U]**
- **Why rank 1:** it is the only venue where (a) an almost exactly analogous paper is in
  print, (b) the format tolerates P3's full evidentiary apparatus — three seeds, bootstrap
  CIs, negative control, matched-capacity PCA arm, placebo deflation — without
  amputation, and (c) the audience already argues about whether interventional data
  improves representations.

### Rank 2 — *Nature Methods* (Brief Communication, or Analysis)

- **Negative results:** **the strongest evidence found anywhere in this brief.** Two
  independent verifications:
  1. **Precedent.** Ahlmann-Eltze et al. 2025 is a Brief Communication (JATS `<subject>`
     tag reads `Brief Communication`) whose title contains "does not yet outperform". **[V]**
  2. **Policy.** The Nature Methods article-types page (fetched from
     `https://www.nature.com/nmeth/content`) defines **Brief Communication** as, verbatim:
     *"A concise report describing potentially groundbreaking yet preliminary method or
     tool developments, highly practical tweaks to an existing method or tool, software
     platforms, resources of broad interest, and **technical critiques of widely used
     methodologies**."* **[V]** That final clause is a standing invitation for P3's genre.
- **Article types and limits (all verbatim from the same page) [V]:**
  - **Brief Communication** — abstract ≤70 words; main text 1,200 words (up to 1,600 at
    editorial discretion) **including abstract, references and figure legends**; max 2
    display items (3 at discretion); ~20 references.
  - **Analysis** — *"A report presenting comprehensive performance comparisons of
    established, related methods or tools, of key importance to a field of research."*
    Abstract ≤150 words; main text 3,000 words (up to 5,000 at discretion); up to 6
    display items; ~50 references.
  - **Article** — *"A report describing a novel method or tool…"* Same limits as Analysis.
- **Open access / APC:** not retrieved. **[U]**
- **Time to first decision:** not retrieved. **[U]**
- **Why rank 2 and not rank 1:** the 1,200-word Brief Communication is brutal for a paper
  whose entire force is in its controls — the exemplar survives it only by pushing ten
  Extended Data figures behind two main figures (§3.1). And *Analysis* is defined around
  comparisons of **"established, related methods or tools"**; P3 compares two supervision
  *targets* inside one architecture, which is a weaker literal fit. **[I]** High reward,
  meaningfully higher rejection risk.

### Rank 3 — *Nature Machine Intelligence*

- **Scope fit:** ML-for-science; comfortable with representation-learning claims. **[I]**
- **Negative results:** three separate routes, all verified from the journal's
  article-types page (`https://www.nature.com/natmachintell/content`) **[V]**:
  - **Analysis** — *"A new analysis of existing data or describes new data obtained in a
    comparative analysis that leads to novel and arresting conclusions."* 3,500 words,
    ≤6 display items, abstract 100–150 words. **This is the single best literal
    article-type match for P3 found in this entire survey.**
  - **Matters Arising** — *"Exceptionally interesting and timely scientific comments and
    clarifications on original research papers published in Nature Machine Intelligence."*
  - **Reusability Reports** — *"Articles that specifically test the robustness and
    reusability of previously published code that supported the findings in papers."*
    A Crossref query on `issn:2522-5839` for the title phrase returns **24 records**, at
    least twelve of which are literally titled `Reusability report: …`, running roughly
    4–9 journal pages each. **[V]** This is a real, regularly-used channel — though it is
    scoped to re-testing *someone else's published code*, which P3 is not. **[I]**
- **Published negative precedent:** Boettcher 2022 (Nat Mach Intell 5:24–25), *"Inability
  of a graph neural network heuristic to outperform greedy algorithms…"* **[V]**; and the
  Boiarsky et al. / reply pair, *"Deeper evaluation of a single-cell foundation model"*
  (Nat Mach Intell 6:1443–1446, 2024, DOI
  [10.1038/s42256-024-00949-w](https://doi.org/10.1038/s42256-024-00949-w)) with *"Reply
  to: …"* (DOI [10.1038/s42256-024-00948-x](https://doi.org/10.1038/s42256-024-00948-x)). **[V]**
- **Open access / APC / decision time:** not retrieved. **[U]**

### Rank 4 — MELBA (*Machine Learning for Biomedical Imaging*)

- **Scope fit:** this is the venue that keeps P3 *inside the computational-pathology
  community* rather than exporting it to comp-bio. **[I]**
- **What was verified** (from `https://www.melba-journal.org/`): *"A web-based journal
  devoted to the free and unrestricted access of high quality articles in the broad field
  that bridges machine learning and biomedical imaging."* The journal charges **no
  publication fees**; there is a **$10 Scholastica submission charge** which the editors
  state they are *"actively working on removing"*. The site's framing is *"you wrote it,
  the community reviewed it, we publish it – no hidden charges and you own your own
  publication."* **[V]**
- **Negative results:** **no explicit policy found** on the pages retrieved. **[U]** It
  publishes benchmarking and evaluation work — e.g. the Ammeling et al. mitotic-figure
  benchmark (vol. 3, 2026) and *"Exploring Fairness and Performance Drivers Across
  State-of-the-Art Pulmonary Nodule Detection Algorithms"* (DOI
  [10.59275/j.melba.2025-6838](https://doi.org/10.59275/j.melba.2025-6838), 2025). **[V]**
- **Length limits / decision time:** not retrieved. **[U]**
- **Why rank 4:** low friction, no APC, right community — but low citation reach, and the
  Ganz→Ammeling episode in §1.1 is a caution that even here the *negatively-framed*
  version was the one that stayed in the workshop.

### Rank 5 — TMLR (*Transactions on Machine Learning Research*)

- **Negative results:** **the most explicitly hospitable acceptance criteria of any venue
  surveyed**, though never using the words "negative results". From
  `https://jmlr.org/tmlr/acceptance-criteria.html`, the two criteria are that *"claims made
  in the submission are supported by accurate, convincing and clear evidence"* and that
  *"some individuals in TMLR's audience [would] be interested in the findings of this
  paper"*. The page explicitly disclaims requirements for state-of-the-art results,
  methodological novelty, or demonstrated significance, stating that *"novelty of the
  studied method is not a necessary criteria for acceptance."* **[V]** The journal
  homepage frames this as prioritising *"technical correctness over subjective
  significance"*. **[V]**
- **Process:** rolling submissions, double-blind, hosted on OpenReview with open review;
  variable manuscript length. TMLR operates an **"Outstanding Certification"** and has
  joined the **NeurIPS/ICML/ICLR Journal-to-Conference Track**, so a TMLR acceptance can be
  presented at a major conference. **[V]**
- **Length limits / APC / stated decision time:** not stated on the pages retrieved. **[U]**
- **Why only rank 5:** the criteria are tailor-made for P3, but the audience is machine
  learning, not pathology or cancer biology. A refuted hypothesis about Perturb-seq
  dictionaries versus Hallmark gene sets accrues most of its value in front of readers who
  care which of those two things is true. **[I]**

### Rank 6 — *PLOS ONE* (floor) / *GigaScience* (floor)

- **PLOS ONE — the only venue surveyed with an explicit, quotable negative-results
  policy.** From `https://journals.plos.org/plosone/s/criteria-for-publication`, verbatim:
  ***"In keeping with our mission to publish all valid research, we consider negative and
  null results."*** **[V]** Its seven criteria require technical soundness, conclusions
  *"supported by the data"*, and data availability — but **no** novelty or significance
  threshold. **[V]** Replication is permitted with *"a sound scientific rationale"*. **[V]**
  APC and decision time not retrieved. **[U]**
- **GigaScience** — from `https://academic.oup.com/gigascience/pages/instructions_to_authors`:
  article types are **Commentary, Data Note, Research, Review, Technical Note, Brief
  Communication**, and — the useful part — ***"Criteria for publication are
  reproducibility, usability and utility, rather than subjective assessment of
  'impact'."*** **[V]** CC BY 4.0. **[V]** No explicit negative-results statement was
  found. **[U]** APC not stated on that page. **[U]**
- **Why floor and not target:** both will take the paper; neither will get it read by the
  people whose hypothesis it refutes.

### Venues to flag

| Venue | Flag | Evidence |
|---|---|---|
| **PLOS Computational Biology** | **Discourages.** Its journal-information page requires *"exceptional significance"*, *"high importance to the field"*, *"originality"* and *"innovation"*. A refutation is none of those on a literal reading. | `https://journals.plos.org/ploscompbiol/s/journal-information` **[V]** — *note the irony: Venet et al. (§3.3), the best null-result exemplar in this document, is a PLOS Comp Biol paper. The 2011 bar was evidently different.* **[I]** |
| **MICCAI** | **Poor fit, not a prohibition.** Theme *"From Algorithms to Clinical Translation"*; solicits *"translational solutions to real-world clinical challenges"*; no tracks, no article types, no mention of negative/validation/benchmarking work. | MICCAI 2026 CFP **[V]** |
| **MICCAI — timing** | **MICCAI 2026 is closed.** Abstracts were due 2026-02-12 and papers 2026-02-26; the conference runs 2026-09-27 → 10-01 in Strasbourg. Next realistic entry is MICCAI 2027 (≈Feb 2027 deadline **[I]**). | MICCAI 2026 CFP **[V]** |
| **Bioinformatics (OUP)** | **No suitable article type.** Original Paper ≤7 pages/~5,000 words; Application Note ≤4 pages/~2,600 words; Reviews 3–8 pages; Letters to the Editor. The guidelines *"do not explicitly identify any article type designated for critical or negative findings"*. Fully open access, **APC $3,625 USD**. | `https://academic.oup.com/bioinformatics/pages/author-guidelines` **[V]** |
| **Nature Communications** | **Precedent weaker than it looks — do not cite it as a negative-results venue without checking.** The one hit retrieved, *"Limitations of representation learning in small molecule property prediction"* (Dias AL, Bustillo L, Rodrigues T; Nat Commun **14**:6394, 2023, DOI [10.1038/s41467-023-41967-3](https://doi.org/10.1038/s41467-023-41967-3)), is **1,319 body words, 15 references, 0 figures** and was published the *same day* as the research paper it discusses (DOI [10.1038/s41467-023-41948-6](https://doi.org/10.1038/s41467-023-41948-6)). It is a companion **Comment**, not a negative research article. | Crossref + PMC10575963 full text **[V]** |

---

## 3. P3 — structural exemplars

Three, all with section order, counts and quotations read directly from full text.

### 3.1 Exemplar A — the one to model on

> **Ahlmann-Eltze C, Huber W, Anders S. "Deep-learning-based gene perturbation effect
> prediction does not yet outperform simple linear baselines." *Nature Methods*
> **22**(8):1657–1661 (2025). DOI [10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6).
> Open access; PMC12328236; PMID 40759747.** **[V]**

**Article type:** Brief Communication (JATS `<subject>` = `Brief Communication`). **[V]**

**Why structurally analogous:** a sophisticated, well-funded class of models is compared
against *deliberately* simple baselines on a held-out prediction task, and loses; the
paper *is* that finding. It shares P3's exact defensive geometry — the authors must
simultaneously prove they ran the sophisticated method competently and that the simple
baseline is not secretly cheating. It also shares the domain furniture: Perturb-seq-style
perturbation data (Norman, Replogle, Adamson), held-out perturbations, linear/PCA-flavoured
comparators.

**Actual section order [V]:** `Main` → `Methods` (`Data`; `Software versions and
parameters`; `Double perturbation benchmark setup`; `Single perturbation benchmark setup`;
`Evaluation metrics`) → `Reporting summary` → `Online content` → `Supplementary
information`.

**Length and display items [V]:** ~3,950 body words total; 34 references; **12 figure
elements** — 2 main figures plus 10 Extended Data figures. Zero tables. The Main section
runs 25 paragraphs.

**How it handles our shared rhetorical problem [all quotes verbatim from PMC12328236]:**

- **How it opens** — with the *claim*, not the refutation. Three sentences of context, then
  the target: *"Two recent models—scGPT and scFoundation—claim to be able to predict gene
  expression changes caused by genetic perturbations."* The word "claim" does all the work.
  The negative finding is not previewed in paragraph 1; the reader is walked to it.
- **Declaring the baseline's simplicity as a feature, up front** — *"we benchmarked the
  performance of these models against GEARS and CPA and against deliberately simplistic
  baselines."* "Deliberately" appears in paragraph 2 and recurs at the close. This is
  precisely the move P3 needs for ordinary PCA at matched capacity.
- **The core asymmetry argument — the single most adaptable sentence in this document:**
  > *"As our deliberately simple baselines are incapable of representing realistic
  > biological complexity, yet were not outperformed by the foundation models, we conclude
  > that the latter's goal of providing a generalizable representation of cellular states
  > and predicting the outcome of not-yet-performed experiments is still elusive."*

  The logic: the baseline's *inadequacy* is what makes losing to it damning. P3's version
  writes itself — a 128-dimensional interventional dictionary that cannot beat ordinary PCA
  at matched capacity has failed at the one thing that was supposed to distinguish it.
- **Foreclosing "you just did it wrong" — by auditing the original papers' own comparisons:**
  > *"The publications that presented GEARS, scGPT and scFoundation included comparisons
  > against GEARS and CPA and against a linear model. Some of these comparisons may have
  > happened to be particularly 'easy'. … The linear model used in scGPT's benchmark
  > appears to have been set up such that it reverts to predicting no change over the
  > control condition for any unseen perturbation."*

  Rather than defending their own competence, they show the prior positive results came
  from a weak comparator. Note the hedging — "may have happened to be", "appears to have
  been" — the accusation is made without ever being asserted.
- **Converging evidence as armour:** *"Our results are in line with previously published
  benchmarks that assessed the performance of foundation models for other tasks and found
  negligible benefits compared to simpler approaches. … Since the release of our paper as
  a preprint, several other benchmarks were released that also show that deep learning
  models struggle to outperform simple baselines."* They also **adopt a critic's better
  baseline mid-review**: *"following the preprints by Kernfeld et al. and Csendes et al.
  that appeared while this paper was in revision"* — visibly conceding on a sub-point to
  buy credibility on the main one.
- **Limitations, stated flatly and early-ish:** *"One limitation of our benchmark is that
  we used only four datasets. We chose these as they were used in the publications
  presenting GEARS, scGPT and scFoundation."* The justification is airtight — the datasets
  were chosen *by the people being critiqued*.
- **The constructive close (never triumphalist):** *"Deep learning is effective in many
  areas of single-cell omics. However, prediction of perturbation effects still remains an
  open challenge, as our present work shows. We expect that increased focus on performance
  metrics and benchmarking will be instrumental to facilitate eventual success…"*
- **The title hedge:** *"does **not yet** outperform"*. "Yet" converts an attack into a
  progress report. Cheap, and it works.

**Caveat worth internalising:** this paper still drew the titled bioRxiv rebuttal in §1.4.
Best-in-class execution did not prevent the counter-attack; it just meant the counter-attack
had to be about metric choice rather than competence. That is the win condition.

### 3.2 Exemplar B — the "most favourable setting" defence

> **Kedzierska KZ, Crawford L, Amini AP, Lu AX. "Zero-shot evaluation reveals limitations
> of single-cell foundation models." *Genome Biology* **26**(1):101 (2025). DOI
> [10.1186/s13059-025-03574-x](https://doi.org/10.1186/s13059-025-03574-x). Open access;
> PMC12007350.** **[V]**

**Why structurally analogous:** the contribution is that an *evaluation choice* — refusing
to fine-tune — exposes a weakness invisible under the standard protocol. P3's contribution
is the mirror image: an *evaluation choice* (40 targets neither arm was trained on, matched
capacity by construction) exposes that the supervision target everyone expected to win,
loses. Both papers argue that the protocol is the contribution.

**Actual section order [V]:** `Background` → `Results and discussion` → `Conclusions` →
`Methods` (`Models and baselines`; `Datasets`; `Evaluation metrics`; `Biological
preservation scores`; `Batch mixing scores`; `Reconstructing gene expression`) →
`Supplementary information`.

Note the merged **`Results and discussion`** and the separate **`Conclusions`** — the
interpretation lives beside the evidence, and `Conclusions` is used for scope, limitations
and forward work. For a paper defending a null, this is a good shape: it prevents a
Discussion section that reads as a long apology.

**Length and display items [V]:** ~4,918 body words; **2 main figures**; 0 tables; 45
references. Everything else — 7+ supplementary figures, 6 supplementary tables — is in
additional files.

**How it handles our shared rhetorical problem [verbatim, PMC12007350]:**

- **Pre-emptive "we gave them every advantage" — steal this construction wholesale:**
  > *"we test this claim and show that even with a set of benchmarks representing **the
  > most favorable setting** where datasets consist of tissues and are generated using
  > technologies similar to those used to pretrain these models …, both Geneformer and
  > scGPT underperform simpler methods."*

  The "you just did it wrong" objection is answered *in the Background*, before any result
  is shown. P3's equivalent: both arms matched by construction on everything but the
  supervision target, evaluated on targets neither arm saw.
- **Naming whose claim is being tested, with a citation:** *"Our evaluations are motivated
  by the authors' claims that their proposed models not only generate robust cell
  embeddings but also exhibit strong capabilities for generalizing to unseen datasets."*
  The paper is not attacking a field; it is testing a specific, attributed proposition.
- **Splitting the failure into competing mechanisms, then testing them** — the direct
  analogue of P3's effective-rank analysis contradicting a capacity explanation:
  > *"We propose two hypotheses as to why Geneformer and scGPT underperform zero-shot
  > compared to the tested baselines. First, it could be that the masked language model
  > pretraining framework … does not produce useful cell embeddings. Second, it could be
  > that scGPT and Geneformer have failed to learn the pretraining task."*

  They then evaluate reconstruction performance to discriminate. A negative result that
  *diagnoses* rather than merely reports is much harder to dismiss.
- **The killer control, understated:** *"scGPT without embeddings underperforms against a
  naive baseline of predicting the mean."*
- **Limitations that hand the reader a fix:** *"One limitation of our analyses was that
  some of our evaluation datasets were used in pretraining, confounding our ability to say
  if performance trends are general…"* — followed immediately by a constructive proposal
  (*"we propose the creation of benchmark tasks and datasets reserved exclusively for model
  evaluation that should never be used to pretrain any future model"*).
- **Explicit generosity to the target:** *"some more recent advances occurring at time of
  our work have already sought to address limitations exposed by our analysis"*, and a
  closing that credits the field — *"foundation models hold significant promise … The
  challenges identified in this study underscore the necessity of a meticulous
  evaluation."*
- **Scope discipline:** *"Recognizing the rapid advancement of the field, our study adopts
  a focused approach."* Three seeds and 40 targets is a focused approach, not a small one.
  Say it that way.

### 3.3 Exemplar C — the random-control paper, and P3's secondary finding

> **Venet D, Dumont JE, Detours V. "Most Random Gene Expression Signatures Are
> Significantly Associated with Breast Cancer Outcome." *PLoS Computational Biology*
> **7**(10):e1002240 (2011). DOI
> [10.1371/journal.pcbi.1002240](https://doi.org/10.1371/journal.pcbi.1002240). Open
> access; PMC3197658.** **[V]**

**Why structurally analogous — and it is the closest match in this document to P3's
*second* claim:** the paper's core result is that **random gene sets reproduce the signal
attributed to curated, biologically-motivated signatures**, and its second move is a
**deflation against a proliferation metagene (meta-PCNA)** showing the shared signal is
largely proliferation. P3 reports (a) most per-target "pathway" signal is reproducible
with random gene sets, and (b) a placebo-controlled deflation shows nothing is special
about proliferation in the arm gap. Venet et al. is the canonical prior work for move (a)
and the methodological ancestor of move (b) — **P3 must cite this paper regardless of which
exemplar it models, and should state explicitly how its placebo-controlled deflation
differs from meta-PCNA adjustment.**

**Actual section order [V]:** `Introduction` → `Results` (five *declarative-sentence*
subheadings, each stating its finding — see below) → `Discussion` → `Methods` (`Software
setup`; `Code and data availability`; `Expression data`; `Literature signatures`;
`Meta-PCNA index`; `Adjusting data for the meta-PCNA index`; `Association of signatures
with outcome`) → `Supporting Information`.

The Results subheadings are themselves arguments — e.g. *"Most signatures not biologically
related to cancer are statistically associated with breast cancer outcome"*, *"Most
published breast cancer signatures are not more strongly associated with breast cancer
outcome than sets of random genes"*, *"Results are reproducible across cohorts and
end-points"*. **Adopt this device.** A reader skimming P3's Results contents should be able
to reconstruct the whole argument, including the negative control and the reproducibility
across seeds.

**Length and display items [V]:** ~6,303 body words; **6 figures**; 1 table; 68 references.

**How it handles our shared rhetorical problem [verbatim, PMC3197658]:**

- **Opens on the *inferential practice*, not the result.** The Introduction spends five
  paragraphs laying out the three-step argument the field uses (characterise a mechanism →
  derive a marker → show the marker correlates with outcome) before touching data. The
  paper is framed as auditing an *argument form*, which is far more durable than auditing a
  method.
- **A lay analogy carries the confounding argument** — memorable and disarming:
  > *"We may find that the number of TV sets per household is positively correlated with
  > longer life expectancy. This, of course, does not imply that TV sets improve health."*

  P3 has an equivalent available: an interventional dictionary can be genuinely causal in
  its source experiment and still be a worse *supervision target* than a curated gene set.
- **Naming the missing control as the field's error, not the authors' cleverness:**
  > *"Few studies using the outcome-association argument present negative controls to check
  > whether their signature of interest is indeed more strongly related to outcome than
  > signatures with no underlying oncological rationale."*
- **Reframing the statistical question — the sentence P3's negative control section should
  be built around:**
  > *"Nominal p-values do not answer the appropriate statistical question: the question is
  > not whether a given set of genes is related to survival, but whether it is **more**
  > related to survival than random sets of genes."*
- **Refusing the over-claim, conspicuously** — the model for P3's placebo-deflation
  paragraph:
  > *"Yet—we cannot stress this enough—we have not shown that proliferation is a core
  > driving force behind breast cancer progression."*

  The em-dashed interjection is doing reviewer-management: it pre-empts the reading that
  would make the paper attackable.
- **Ring-fencing the damage:** *"Our study questions the biological interpretation of the
  prognostic value of published breast cancer signatures, but has no bearing on their
  usefulness in the clinic: a marker may be accurate without yielding interesting
  biological insight regarding the mechanism of disease progression."* P3's analogue:
  refuting the *supervision-target* hypothesis says nothing about whether Perturb-seq is
  useful for the purposes it was designed for.
- **Enumerated conclusion.** The Discussion ends *"In conclusion, we have shown that 1)…
  2)… 3)… 4)…"* — four numbered claims, each independently defensible. For a paper that
  will be attacked, a numbered conclusion forces critics to engage a specific item.

---

## 4. P3 — the single exemplar to model most closely

**Model P3 on Ahlmann-Eltze, Huber & Anders (§3.1), *Nature Methods* 2025.** It is the
closest available match on all four axes that matter: the *contribution type* (a
predeclared head-to-head in which the sophisticated arm loses to a deliberately simple one,
and the paper is about that), the *defensive geometry* (the authors must prove both that
they ran the sophisticated arm competently and that the simple arm is not cheating — P3's
exact burden with the interventional dictionary versus matched-capacity PCA), the *domain
furniture* (Perturb-seq perturbation data, held-out perturbations, linear comparators), and
the *outcome* (it was published in a top methods journal, with "not" in the title, and
survived a public titled rebuttal). Its central rhetorical device — *"As our deliberately
simple baselines are incapable of representing realistic biological complexity, yet were
not outperformed by the foundation models, we conclude…"* — is the argument P3 needs to
make about ordinary PCA, and its Extended-Data strategy (2 main figures carrying the
headline, 10 Extended Data figures carrying the controls, robustness and cost analyses) is
exactly how P3 fits three seeds, bootstrap CIs, a random-gene-set negative control, an
effective-rank capacity analysis and a placebo-controlled deflation into a short-format
paper without losing any of them. **Borrow two things from the other exemplars:** Kedzierska
et al.'s "most favourable setting" sentence, moved into P3's Background so the "you just did
it wrong" objection is answered before the first result; and Venet et al.'s
declarative-sentence Results subheadings, plus its numbered conclusion, so a skimming
reviewer cannot miss that the negative control and the capacity control were run. Cite
Venet et al. explicitly as the ancestor of P3's random-gene-set finding — claiming that
result as novel is the fastest route to a hostile review.

---

## 5. P4 — venue classes

P4 is early (one of five gate conditions met), so this section gives **venue classes with
verified precedent**, not a submission plan. Each class is evidenced by papers retrieved
and independently re-verified against Crossref during compilation.

### The finding that should shape P4's positioning

> **No retrieved paper ships a promptable natural-language interface over a biological
> atlas whose headline feature is *refusing to answer* uncertified queries.** **[V — as a
> negative over the searches run; absence of evidence, not proof of absence]**

The two halves exist separately and both are publishable at the top of the field. They
have not been joined. Concretely:

- The **interface half** is proven. *CZ CELLxGENE Discover* is a 15-page NAR Database-issue
  paper whose entire contribution is a queryable platform. But a full-text search of its
  body for `trust` and `caveat` returned **zero hits** — the resource genre currently makes
  **no epistemic warranty at all**. **[V]**
- The **abstention half** is proven, in clinical AI rather than in resources: CoDoC at
  *Nature Medicine*, popV at *Nature Genetics*, Leibig et al. at *Scientific Reports*.
- Only **popV** does both in one paper — and it is a method that emits a *score*, not an
  interface that operates a *gate*. **[V]**

That gap is P4's claim to novelty and also its principal risk: a genre that has not formed
may not be a gap so much as an absence of demand. Treat the white space as a hypothesis to
be tested at the gate conditions, not as a result.

### The foil P4 must differentiate against — and it is very recent

> **Schaefer M, Peneder P, Malzl D, Lombardo SD, Peycheva M, Burton J, Hakobyan A, Sharma V,
> Krausgruber T, Sin C, Menche J, Tomazou EM, Bock C. "Multimodal learning enables
> chat-based exploration of single-cell data." *Nature Biotechnology*, published online
> 2025-11-11. DOI [10.1038/s41587-025-02857-9](https://doi.org/10.1038/s41587-025-02857-9).
> PMID 41219484. Not open access; no PMCID, so no section structure available.** **[V —
> re-verified via Crossref: Nature Biotechnology, 13 authors, first author Moritz Schaefer]**

This is, verbatim from its abstract, a system where *"This embedding informs a large
language model that answers user-provided questions about cells and genes in
natural-language chats"*, and where the authors *"integrate a CellWhisperer chat box with
the CELLxGENE browser, allowing users to interactively explore gene expression through a
combined graphical and chat interface."* **[V]** Its abstract contains **no** uncertainty,
abstention, calibration, or refusal language; its only reliability framing is a
zero-shot-prediction benchmark. **[V]**

**This is exactly the artefact P4's governing principle says is worse than nothing** — a
fluent chat interface over an uncertified atlas — and it is in *Nature Biotechnology*.
P4 cannot be positioned as "the first chat interface over an atlas". It must be positioned
as **the certificate layer that the existing chat interfaces lack**, with CellWhisperer as
the named, dated, high-profile instance of the problem.

### The classes

| Class | Venues | Verified precedent | What it proves for P4 |
|---|---|---|---|
| **1. Nature-family methods / biotech / genetics** — *strongest overall* | Nature Methods, Nature Biotechnology, Nature Genetics | HTAN DCC, Nat Methods **22**:664–671 (2025), DOI [10.1038/s41592-025-02643-0](https://doi.org/10.1038/s41592-025-02643-0); Vitessce, Nat Methods **22**:63–67 (2025), DOI [10.1038/s41592-024-02436-x](https://doi.org/10.1038/s41592-024-02436-x); CellWhisperer (above); popV, Nat Genet **56**:2731–2738 (2024), DOI [10.1038/s41588-024-01993-3](https://doi.org/10.1038/s41588-024-01993-3) | Publishes **both** halves, and in popV both in one paper. Vitessce shows the short format: ~4,400 body words, 2 figures, 54 refs. **[V]** |
| **2. Database / resource issues** | Nucleic Acids Research; GigaByte | CZ CELLxGENE Discover, NAR **53**:D886–D900 (2025 Database issue; online 2024-11-28), DOI [10.1093/nar/gkae1142](https://doi.org/10.1093/nar/gkae1142); *Portable-CELLxGENE*, GigaByte (2025), DOI [10.46471/gigabyte.151](https://doi.org/10.46471/gigabyte.151) | Strongest precedent for **interface-as-contribution**; weakest for trustworthiness (zero `trust` hits in CELLxGENE body). **[V]** The class where P4's certificate layer is genuinely new. |
| **3. Bioinformatics applications / brief-report track** | Bioinformatics (OUP) | UCSC Cell Browser, DOI [10.1093/bioinformatics/btab503](https://doi.org/10.1093/bioinformatics/btab503) (`pubType: brief-report`); AI-HOPE conversational agent, DOI [10.1093/bioinformatics/btaf359](https://doi.org/10.1093/bioinformatics/btaf359) (`pubType: brief-report`) | An LLM natural-language analysis agent is publishable here as a brief report. Note AI-HOPE's only trust claim is *"its closed-system design prevents clinical data leakage"* — security, not epistemic certification. **[V]** Fallback if P4 lands as a tool, not a resource. |
| **4. Clinical AI / digital medicine** — *the only class where refusal is a title-level claim* | Nature Medicine; npj Digital Medicine; Scientific Reports; JAMIA Open; Artificial Intelligence in Medicine; Statistical Methods in Medical Research; PSB | CoDoC, Nat Med **29**:1814–1820 (2023), DOI [10.1038/s41591-023-02437-x](https://doi.org/10.1038/s41591-023-02437-x); *When silence is safer*, npj Digit Med (2026-06-16), DOI [10.1038/s41746-026-02882-1](https://doi.org/10.1038/s41746-026-02882-1); Leibig et al., Sci Rep **7**:17816 (2017), DOI [10.1038/s41598-017-17876-z](https://doi.org/10.1038/s41598-017-17876-z) | Abstention framed as a capability, not a caveat, at the top clinical venue. **[V]** |
| **5. Cell Press consortium flagship** | Cell | Rozenblatt-Rosen O, Regev A, Oberdoerffer P, Nawy T, Hupalowska A, Rood JE, et al. *The Human Tumor Atlas Network: Charting Tumor Transitions across Space and Time at Single-Cell Resolution.* Cell (2020), DOI [10.1016/j.cell.2020.03.053](https://doi.org/10.1016/j.cell.2020.03.053) | Consortium-launch format. **[I]** P4 does not fit this unless positioned as a collective effort at HTAN scale. |

**Class 4 supplies P4's governing-principle citation.** *When silence is safer* states, in a
Nature-portfolio venue, that *"Large language models (LLMs) are designed to generate answers
to user prompts, which often drives them to respond even when uncertainty is high,
information is incomplete, or a refusal would be more appropriate"*, and that *"confidently
stated but inaccurate medical advice can cause significant harm, making the ability to
abstain especially important."* **[V]** That is "you cannot prompt what you cannot certify",
already in print. **Cite it — it converts P4's principle from a manifesto into an
instantiation of an established design norm.**

---

## 6. P4 — structural exemplars

### 6.1 Exemplar A — the one to model on

> **Ergen C, Xing G, Xu C, Kim M, Jayasuriya M, McGeever E, Oliveira Pisco A, Streets A,
> Yosef N. "Consensus prediction of cell type labels in single-cell data with popV."
> *Nature Genetics* **56**:2731–2738 (2024). DOI
> [10.1038/s41588-024-01993-3](https://doi.org/10.1038/s41588-024-01993-3). PMID 39567746;
> PMC11631762.** **[V — re-verified via Crossref: 9 authors, first author Can Ergen]**

**Why structurally analogous:** popV returns a label **plus a score**, and the paper's
headline value is not accuracy but that low-score regions are **surfaced for human scrutiny
rather than answered confidently**. That is P4's certified/uncertified split, one axis at a
time, already validated at a Nature-family journal. Critically, popV shows how to make
"we decline on part of the output" read as the *product* rather than as a shortfall.

**Actual section order [V]:** `Main` → `Results` (`Overview of popV`; `PopV prediction score
discriminates high- and low-quality annotations`; `PopV provides useful label transfer in
case of drastic differences in cellular composition`) → `Discussion` → `Methods` (28
sub-headings, including `Consensus voting`, `Evaluation metrics`, `Ablation experiment`,
`Statistics and reproducibility`, `Reporting summary`).

**Note the second Results heading.** It is a declarative sentence asserting that *the score
separates trustworthy from untrustworthy output*. **P4's equivalent heading is the whole
paper**: "the certificate discriminates answerable from unanswerable axes". Make it a
Results heading, not a Discussion aspiration.

**Length and display items [V]:** ~9,069 body words; **3 figures**; 0 tables; 35
references. The figure count is strikingly low for a Nature-family paper — the load sits in
supplementary. P4 should plan the same way.

**Verbatim phrases to adapt [V]:**
- ABSTRACT: *"Existing methods for transferring cell-type labels lack proper uncertainty
  estimation for the resulting annotations, limiting interpretability and usefulness."*
  — the gap statement. P4's version names fluent-but-uncertified interfaces.
- ABSTRACT: *"popV confidently annotates the majority of cells while highlighting cell
  populations that are challenging to annotate by label transfer."*
  — **the single most adaptable sentence for P4.** Note the shape: *confident on most,
  explicitly flagged on the rest*. It never uses the words "fails", "cannot" or "limitation".
- BODY: *"it is essential that annotation methods highlight areas of uncertainty that
  require expert knowledge input."*
- BODY: *"it is crucial for automatic cell-type annotation pipelines to highlight areas of
  uncertainty that may require manual scrutiny, balance the specificity of predictions with
  accuracy and be easily accessible and usable."*
  — the three-part design requirement. P4 can state its gate conditions in this register.
- BODY: *"We emphasize that intrinsic and extrinsic uncertainty are two complementary
  measurements essential to quantifying the performance of a set of cell annotation tools."*

### 6.2 Exemplar B — the resource-paper skeleton

> **CZI Cell Science Program; Abdulla S, Aevermann B, Assis P, Badajoz S, Bell SM, Bezzi E,
> Cakir B, Chaffer J, Chambers S, Cherry JM, et al. (54 authors). "CZ CELLxGENE Discover: a
> single-cell data platform for scalable exploration, analysis and modeling of aggregated
> data." *Nucleic Acids Research* **53**:D886–D900 (2025 Database issue; published online
> 2024-11-28). DOI [10.1093/nar/gkae1142](https://doi.org/10.1093/nar/gkae1142). PMID
> 39607691; PMC11701654.** **[V — re-verified via Crossref: 54 authors, consortium first
> author]**

**Why structurally analogous:** the canonical "we shipped a queryable interface over a very
large biological representation, and the paper *is* the interface's warrant". It proves a
platform description is publishable as primary research, and it models the **multi-surface
argument** P4 needs — one representation, several access modes, each justified separately.

**Actual section order [V]:** `Introduction` → `Results` (10 sub-headings, including
`Navigating the CZ CELLxGENE data corpus`; `Scalable tools allow biologists to explore,
query and analyze CZ CELLxGENE data`; `Explorer allows interactive exploration…`; `Gene
Expression allows gene expression queries across the corpus of data`; `Census provides
efficient programmatic access…`) → `Discussion` → `Materials and methods` (13 sub-headings)
→ back matter.

**Structural lesson:** each access surface gets **its own Results sub-heading**, named for
what the user can do. P4 should give the certified query path and the uncertified-axis
display path their own headings — so that the refusal surface is architecturally equal to
the answer surface, not a subsection of it.

**Length and display items [V]:** ~11,375 body words; **6 figures**; 1 table; 72 references.

**Verbatim phrases [V]:**
- BODY: *"CZ CELLxGENE provides an interoperable and dynamic community resource that
  supports a diversity of biological and computational applications."*
- BODY: *"Diverse use cases are enabled by a uniquely large, multi-organ and consistently
  curated data resource."*
- BODY: *"only an estimated 25% of publicly available datasets providing the cell-level
  metadata needed for reuse"* — how a resource paper quantifies the deficiency it fixes.
  P4 needs the analogous number: what fraction of queried axes currently carry no
  certificate.

**How it handles our shared rhetorical problem — by not handling it.** Grep of the full
body for `trust` and `caveat`: **zero hits.** **[V]** The genre's current norm is to
describe capability and stay silent on warranty. P4's contribution is legible precisely
against that silence — quote this absence rather than asserting that resources are
overconfident.

### 6.3 Exemplar C — making the abstention curve the primary evaluation axis

> **Leibig C, Allken V, Ayhan MS, Berens P, Wahl S. "Leveraging uncertainty information from
> deep neural networks for disease detection." *Scientific Reports* **7**:17816 (2017). DOI
> [10.1038/s41598-017-17876-z](https://doi.org/10.1038/s41598-017-17876-z). PMID 29259224;
> PMC5736701.** **[V — re-verified via Crossref: 5 authors, first author Christian Leibig]**

**Why structurally analogous:** the paper's evaluation is organised around **how performance
improves as the system is allowed to decline more cases**. That is the measurement P4 needs
and currently lacks a template for: not "is the answer right", but "what does gating buy
you, as a function of how much you gate". It also solves P4's hardest presentational
problem — how to show that refusing is *valuable* rather than merely *safe*.

**Actual section order [V]:** `Introduction` → `Results` (`Uncertainty rank orders
prediction performance`; `Performance improvement via uncertainty-informed decision
referral`; `Performance improvement for different costs, networks, tasks and datasets`;
`Comparison with alternative uncertainty measures`; `What causes uncertainty?`; `Uncertainty
about unfamiliar data samples`) → `Discussion` → `Methods` → back matter.

Again: **declarative Results headings**, each asserting a finding. And note the sequence —
*it works* → *it works across conditions* → *it beats alternatives* → *here is why it
works*. P4's certificate section can follow that arc exactly.

**Length and display items [V]:** ~10,883 body words; **8 figures**; 2 tables; 75 references.

**Verbatim phrases to adapt [V]:**
- ABSTRACT: *"a physician knows whether she is uncertain about a case and will consult more
  experienced colleagues if needed"* — the human-competence analogy that makes abstention
  read as expertise rather than weakness. P4 should open on something of this shape.
- ABSTRACT: *"we show that uncertainty informed decision referral can improve diagnostic
  performance"* — note **"improve"**. Refusal framed as a performance gain.
- BODY: *"The computed measure of uncertainty allowed us to refer a subset of difficult
  cases for further inspection, resulting in substantial improvements in detection
  performance in the remaining data."*
- BODY: *"The decision referral scenario can serve as a minimal benchmark for comparing
  uncertainty methods."* — **directly reusable**: precedent for making the abstention curve
  the primary evaluation axis rather than a supplementary robustness check.

### 6.4 Secondary references worth citing (verified, not modelled on)

- **CoDoC** — Dvijotham K, Winkens J, Barsbey M, Ghaisas S, Stanforth R, Pawlowski N, et al.
  (30 authors). *Enhancing the reliability and accuracy of AI-enabled diagnosis via
  complementarity-driven deferral to clinicians.* Nat Med **29**:1814–1820 (2023). DOI
  [10.1038/s41591-023-02437-x](https://doi.org/10.1038/s41591-023-02437-x). **[V]**
  No PMCID → no structure available. **[V]** Abstract quotes: *"a system that can learn to
  decide between the opinion of a predictive AI model and a clinical workflow"*; *"reduced
  false positives by 25% at the same false-negative rate, while achieving a 66% reduction
  in clinician workload"* — the template for quantifying what deferral buys.
- **Vitessce** — Keller MS, Gold I, McCallum C, Manz T, Kharchenko PV, Gehlenborg N. Nat
  Methods **22**:63–67 (2025). DOI
  [10.1038/s41592-024-02436-x](https://doi.org/10.1038/s41592-024-02436-x). **[V]**
  Structure: `Main` → `Methods` (20 sub-headings). ~4,405 body words; 2 figures; 54 refs.
  **[V]** Use as the **length/format calibration point** for a short-format Nature Methods
  interface submission.
- **HTAN DCC** — de Bruijn I, Nikolov M, Lau C, Clayton A, Gibbs DL, Mitraka E, et al. (27
  authors). Nat Methods **22**:664–671 (2025). DOI
  [10.1038/s41592-025-02643-0](https://doi.org/10.1038/s41592-025-02643-0). **[V]** Full
  text returned empty from Europe PMC, so **no structure, word count or figure count is
  reported here.** **[U]** Abstract quote **[V]**: *"HTAN data can be accessed through the
  HTAN Portal, explored in visualization tools-including CellxGene, Minerva and
  cBioPortal-and analyzed in the cloud through the NCI Cancer Research Data Commons."* —
  the multiscale-tumour-atlas substrate P4 sits on.
- **Bio-SODA UX** — Sima AC, Mendes de Farias T, Anisimova M, Dessimoz C, Robinson-Rechavi
  M, Zbinden E, Stockinger K. *Bio-SODA UX: enabling natural language question answering
  over knowledge graphs with user disambiguation.* Distributed and Parallel Databases
  **40**:409–440 (2022). DOI
  [10.1007/s10619-022-07414-w](https://doi.org/10.1007/s10619-022-07414-w); PMC9458692.
  **[V]** ~12,595 body words; 8 figures; 6 tables; 48 refs. **[V]** Body quote: *"Since
  natural language questions can be very ambiguous and two users might mean different
  things when asking the same question, it is important that the system helps the user in
  exploring the data and in finding the correct answer to a question."* **[V]**
  **Caveat:** its mechanism is *user disambiguation*, not a statistical certificate — the
  interaction-pattern analogue, not the guarantee analogue. Venue is a database journal.

---

## 7. P4 — the single exemplar to model most closely

**Model P4 on Ergen et al., popV (§6.1), *Nature Genetics* 2024.** It is the only retrieved
paper that carries both of P4's halves in one artefact: a shipped, usable resource *and* a
principled account of what it declines to assert, with the declining presented as the
product rather than as a limitation. Its abstract sentence — *"popV confidently annotates
the majority of cells while highlighting cell populations that are challenging to annotate
by label transfer"* — is the exact register P4 needs for "answers certified axes, displays
uncertified ones as uncertified": confident on most, explicitly flagged on the rest, with no
vocabulary of failure anywhere in the sentence. Its Results architecture is equally
transferable: a declarative heading asserting that the score *discriminates* trustworthy
from untrustworthy output, which for P4 becomes the load-bearing claim that the certificate
discriminates answerable from unanswerable axes — a Results finding to be demonstrated, not
a Discussion aspiration. Its proportions (~9,000 body words, only 3 main figures, heavy
supplementary) show how to keep a resource paper's evidentiary burden without a wall of
display items. **Borrow from the other two:** CELLxGENE's device of giving each access
surface its own named Results sub-heading — so P4's uncertified-axis display gets equal
architectural weight to its answer path, rather than reading as an appendix to it; and
Leibig et al.'s decision-referral curve as the *primary* evaluation axis, since it is the
only retrieved precedent for demonstrating that gating **improves** rather than merely
protects. Finally, note that P4's differentiation problem is now concrete and dated:
CellWhisperer (*Nature Biotechnology*, November 2025) already ships a chat box over
CELLxGENE with no uncertainty language in its abstract at all. P4 is not the first
promptable atlas; it is the certificate layer that the first ones did not have, and every
draft should be written from that position.

---

## 8. Verification ledger — what stayed unverified

Check these before relying on them. Nothing below was guessed at in the body of this
document; each is flagged **[U]** at point of use.

**P3**
1. MICCAI page/length limit — the 2026 and 2025 author-guideline pages both 302-redirect to
   the bare conference root; not retrieved. MICCAI 2026 dates and deadlines **were**
   verified from the CFP.
2. *Medical Image Analysis* guide-for-authors — not fetched. No negative-results article
   type was found via Crossref, but the policy page itself is unread.
3. *Genome Biology* article types, length limits and APC — submission-guidelines page sits
   behind a two-hop Springer identity redirect; not retrieved. Its **precedent** (Kedzierska
   et al.) and open-access status **were** verified.
4. Time to first decision — **not verified for any venue in this document.** No publisher
   page retrieved stated one. Every timing statement you may want is currently missing.
5. APCs — verified only for **Bioinformatics ($3,625 USD)** and **MELBA (no fee; $10
   Scholastica submission charge)**. Nature Methods, Nature Machine Intelligence, Genome
   Biology, PLOS ONE, GigaScience, TMLR: not retrieved.
6. TMLR length limits and review timeline — the acceptance-criteria page gives criteria but
   no timings; the homepage says "shortened review period" without a number.
7. Boettcher 2022 article type — inferred as Matters Arising from its 2-page extent
   (Nat Mach Intell 5:24–25); the article page was not fetched.
8. Raghu et al. *Transfusion* — venue verified only from the arXiv `<arxiv:comment>` field
   reading "NeurIPS 2019". Section order, length and figure count **not** extracted; it is
   cited in §1.2 as precedent only, and is deliberately **not** used as a structural
   exemplar.

**P4**
9. HTAN DCC (PMC12125965), HTAN *Cell* 2020 (PMC7376497) and ClinGen (PMC11984750) all
   returned **empty** full-text XML from Europe PMC. No section structure, word count or
   figure count is reported for any of them.
10. CellWhisperer and CoDoC have **no PMCID** and are not open access — abstracts were
    verified, section structures could not be.
11. NAR Database-issue classification for CELLxGENE — inferred from the `D886–D900`
    pagination, not from a stated article type.
12. The "no paper joins promptable interface + refusal-as-headline" claim in §5 is a
    negative over the searches actually run. It is a strong hint, not a proof.
13. ClinGen was flagged by the discovery pass as an untested lead — a resource whose formal
    output includes an explicit "uncertain" verdict. Worth one hour of follow-up; **not**
    verified here and **not** cited as an exemplar.
