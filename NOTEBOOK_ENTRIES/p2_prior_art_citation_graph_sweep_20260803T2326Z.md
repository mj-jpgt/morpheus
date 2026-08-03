# P2 prior art: citation-graph sweep, closing the `[SEARCH INCOMPLETE]` flag on draft §2.2

UTC 2026-08-03T23:26Z. Branch `research/rebase-vision`.

Task 1 of the P2 assignment. Draft §2.2 carries
`[SEARCH INCOMPLETE — must be closed before submission.]` with the note that the sweep behind it ran
"through the arXiv Atom API and OpenAlex only" and that "a Semantic-Scholar citation-graph sweep of
RankMe's citing papers has not been run". That sweep has now been run. **This entry does not edit
`paper/P2_RANK_DRAFT.md`** — another agent owns it.

---

## 0. Verdict, bad news first

**The paper does not collapse to a replication, but §2.2's central sentence is now false as written
and must be rewritten.** That sentence reads:

> "Every one of these tests rank *across* methods, *across* checkpoints, or in domains other than
> joint-embedding representation selection. None tests the within-method, matched-arm,
> fixed-architecture regime that RankMe reserves for itself…"

Three published papers now falsify it, one of them peer-reviewed and one of them a partial
pre-emption of the **necessity** claim that Task 3 is built on:

1. **Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko (ICASSP 2025)** report rank failing as a
   *selection rule* within a single encoder, and explicitly report **lower-ranked layers
   outperforming higher-ranked ones** — i.e. a published low-rank/high-information instance.
2. **Otero, Mateus & Balestriero (2024)** — already in §2.2, but it is a *selection-rule* failure
   ("cross-validation without labels infeasible"), which §2.2 currently under-states.
3. **Zhang, Jiang, Gao, Willett & Maire (2024)** report a matched-arm anti-correlation:
   "correlation between **low** effective feature rank and downstream task performance".

And on the other side, **three strong defences exist that the draft does not cite at all**, including
one peer-reviewed and **in our own domain** (BMC Genomics, RankMe on DNA language models, concluding
it works) and one running **precisely the hyperparameter-sweep design we claim breaks** and reporting
a power law.

What survives, and is now much better evidenced: across **453 de-duplicated citing works**, there is
**no** within-method matched-arm evaluation of rank as a selection rule in computational pathology or
transcriptomics, and no paper that isolates rank from information content as the explanatory claim.
That is where the contribution should be positioned.

---

## 1. Method and completeness

Semantic Scholar Graph API citation pulls, `/paper/CorpusID:<id>/citations`, plus OpenAlex, Crossref,
arXiv Atom and DBLP for field verification.

| seed | verified record | citers pulled |
|---|---|---|
| **RankMe** | S2 `corpusId 252735059`; externalIds `{ArXiv: 2210.02885, DBLP: conf/icml/GarridoBNL23}`; venue "International Conference on Machine Learning" | **159** |
| **Roy & Vetterli 2007** | S2 `corpusId 12184201`; `{DBLP: conf/eusipco/RoyV07, DOI: 10.5281/ZENODO.40328}` | **581**, of which **325** survived a deep-learning-context filter |
| **LiDAR** | S2 `corpusId 266052843`; `{ArXiv: 2312.04000, DBLP: conf/iclr/Thilak0SDGNSL24}` | **33** |
| **α-ReQ** (added) | S2 `corpusId 286313323`; DOI `10.52202/068431-1281` | **15** |

**Union, de-duplicated by normalised title: 453 unique works.** ~70 abstracts read in full.

### Completeness limits — stated, not papered over

1. **OpenAlex is unusable as a citation source here.** Its two RankMe records return **12 unique
   citers combined**; Semantic Scholar returns **159** for the same paper. A search for a
   PMLR/ICML-2023 RankMe record found none (`title.search:RankMe` → 5 hits, only the 2 arXiv records
   plus unrelated "RankME"/"Rankmed" papers; `raw_author_name.search:Quentin Garrido` +
   `publication_year:2023` → 4 hits, no RankMe). PMLR does not mint Crossref DOIs, which is the
   likely cause. **The earlier sweep behind §2.2 was OpenAlex-based and therefore saw roughly 8% of
   RankMe's citation graph.** That is the root cause of the `[SEARCH INCOMPLETE]` flag and it is now
   fixed.
2. **S2 `/paper/search` was unusable** (HTTP 429 through ~12 backoff attempts). Worked around by
   pulling RankMe's *reference* list to recover Roy & Vetterli's corpus ID. All `/citations` calls
   succeeded with backoff.
3. **Triage is abstract-level only.** No full texts were read. **12 of 159** RankMe citers had no
   abstract in the S2 response and were triaged on title alone. Since LiDAR was nearly missed
   precisely this way, the "noise" bucket is soft and should be treated as a residual risk.
4. **α-ReQ leg is partial** — no arXiv record for it was found (two Atom queries returned
   `totalResults=0`), so a more-cited preprint record may exist unswept.

---

## 2. PRE-EMPTS — must be engaged head-on

### A1. The most dangerous one. **Lead §2.2 with this, do not let a referee find it first.**

**Aldeneh, Thilak, Higuchi, Theobald & Likhomanenko, "Towards Automatic Assessment of Self-Supervised
Speech Models using Rank", ICASSP 2025.**

- DOI **`10.1109/ICASSP49660.2025.10889651`**, arXiv **2409.10787** (v2).
- **Fields I verified myself, from source, independently of the sweep agent:**
  - Crossref `https://api.crossref.org/works/10.1109/ICASSP49660.2025.10889651` →
    title "Towards Automatic Assessment of Self-Supervised Speech Models using Rank"; authors
    `['Zakaria Aldeneh','Vimal Thilak','Takuya Higuchi','Barry-John Theobald','Tatiana Likhomanenko']`;
    container-title "ICASSP 2025 - 2025 IEEE International Conference on Acoustics, Speech and Signal
    Processing (ICASSP)"; published **2025-04-06**.
  - arXiv Atom `https://export.arxiv.org/api/query?id_list=2409.10787` → `arxiv.org/abs/2409.10787v2`,
    published 2024-09-16, same five authors, abstract as quoted below.
- **Verbatim abstract, the two load-bearing sentences:**
  > "The findings indicate rank correlates with downstream performance **within encoder layers**
  > across various downstream tasks and for in- and out-of-domain scenarios. **However, rank does not
  > reliably predict the best-performing layer for specific downstream tasks, as lower-ranked layers
  > can outperform higher-ranked ones.**"

**Why this is the most dangerous item in the sweep — it hits *both* of our contributions:**

- It is a **within-method, matched-arm, selection-rule failure**: one encoder, layers as the matched
  arms, rank picks the wrong one. §2.2's "None tests the within-method … regime" is false.
- **"lower-ranked layers can outperform higher-ranked ones" is a published low-rank/high-information
  instance — the exact configuration Task 3 set out to find.** Our necessity argument is therefore
  *partially pre-empted*, and §5.3 of the companion metrics entry must be reframed as
  *corroborating* Aldeneh et al. in a new regime, not as discovering the configuration.
- **Vimal Thilak is LiDAR's first author**, so this is a negative result about rank published by the
  people who proposed the leading replacement — the same hostile-witness value as Otero/Balestriero.
- It follows the "negative advertised as positive" pattern the assignment warned about: the title and
  the closing sentence sell rank as "a valuable tool for monitoring training progress".

**What is still ours after A1.** Aldeneh et al. vary **layer depth within one trained encoder**. Layers
of a single network are not independent training runs; RankMe's reserved scope is *"different runs of
a given method"*, and a layer sweep is not that. Our matched arms are **separately trained runs**
differing in one objective term, with three seeds each and a reproducibility floor measured on the
statistic. Also: they report rank *correlating* with performance within layers and failing only at
argmax; we report a variance decomposition in which the arm term is not significant at all. The
distinction is real but it is **narrow**, and the paper must state it rather than imply novelty.

### A2. Otero, Mateus & Balestriero (2024) — already cited, but under-stated

arXiv **2410.04289**, v1 2024-10-05, authors `Daniel Otero, Rafael Mateus, Randall Balestriero`
(arXiv Atom verified). No Crossref record for the arXiv DOI.

> "we emphasize the need for better label-free assessments of SSL representations, as **current
> methods like RankMe fail to adequately evaluate representation quality, making cross-validation
> without labels infeasible**."

§2.2 currently files this under "Further negatives — PARTIALLY VERIFIED (abstract only)" and treats
it as out of scope. "Cross-validation without labels infeasible" is a **named selection-rule
failure**, on 250 experiments in which backbone and SSL framework are held fixed while class
imbalance is varied. It is closer to in-scope than the draft allows.

### A3. Zhang, Jiang, Gao, Willett & Maire, "Residual Connections Harm Generative Representation Learning" (2024)

arXiv **2404.10947** (v5), published 2024-04-16, authors verified from arXiv Atom by me.

> "Analyzing the representations learned by our modified residual networks, we find **correlation
> between low effective feature rank and downstream task performance**."

Matched-arm by construction (same MAE / ViT-B/16, one architectural scalar varied), and the *better*
arm — kNN 27.4% → **63.9%**, linear probe 67.8% → 72.7% — is reported as the *lower*-rank one.
**`[ABSTRACT-LEVEL ONLY — the sign must be confirmed in the body before this is cited as a necessity
violation.]`** If the body confirms it, this is a second published low-rank/high-information result
and Task 3's framing weakens further.

### A4. LiDAR itself

Already in §2.2. Noting only that its verbatim *"significantly surpasses naive rank based approaches
in its predictive power of optimal **hyperparameters**"* is within-method matched-arm selection, so
§2.2's claim that none of the prior negatives is in-scope was already strained by a paper the draft
cites.

---

## 3. STRONGEST DEFENCES OF RANK — argue against these, not a straw man

The assignment asked specifically for the best form of the opposing case. **The draft currently cites
none of these three.**

### C1. The design we claim breaks, reported working, with a power law

**Deng, Sun, Dou & Xu, "Unify Variables in Neural Scaling Laws for General Audio Representations via
Embedding Effective Rank", arXiv:2510.10948.**

- **Verified by me from arXiv Atom** (`https://export.arxiv.org/api/query?id_list=2510.10948`):
  `arxiv.org/abs/2510.10948v1`, published **2025-10-13**, authors
  `['Xuyao Deng','Yanjie Sun','Yong Dou','Kele Xu']`. `[UNVERIFIED: peer-reviewed venue — none found;
  cite as an arXiv preprint.]`
- **Verbatim:**
  > "we present a systematic study of scaling laws for general audio representations by utilizing
  > **embedding effective rank (RankMe) as a unifying metric** … allowing us to examine scaling
  > behaviors across a wide hyper-parameter space, including **model size, training data volume,
  > computational budget, architectural configurations**, etc. Our empirical findings reveal a
  > **consistent power-law relationship between RankMe and representation quality**, suggesting that
  > embedding effective rank serves as a **reliable proxy** for assessing and predicting model
  > performance."

This is the strongest counter-case in the corpus: within-method, hyperparameter-swept, matched — the
exact regime we claim breaks — and it reports not merely a correlation but a *power law*. **The paper
must address it directly.** The available lines of attack, in order of strength: (i) it sweeps
*capacity-like* variables (model size, data volume, compute) which move rank and quality together for
reasons that have nothing to do with rank being a quality proxy — our arms are matched **on capacity**
and differ only in an objective term, which is the case that separates the two; (ii) it reports no
seed-level reproducibility floor for RankMe, and our §4 measurement is that the seed term dominates;
(iii) audio, not pathology.

### C2. A peer-reviewed RankMe-works result **inside our own domain**

**Awasthi, Mend Mend Arachchige & Zhu, "Unsupervised evaluation of pre-trained DNA language model
embeddings", BMC Genomics (2025).**

- **Verified by me from Crossref** (`https://api.crossref.org/works/10.1186/s12864-025-11913-2`):
  title as given; authors `['Raghav Awasthi','Gayan Samuditha Mend Mend Arachchige','Xiaofeng Zhu']`;
  journal **BMC Genomics**; published **2025-08-01**; **volume 26, issue 1, article number 710**.
  PMID `40751178`, PMC `12315385` (reported by the sweep; I verified the Crossref fields only).
- **Verbatim:**
  > "We propose a framework to evaluate DLM embeddings using unsupervised numerical linear
  > algebra-based metrics **RankMe, NESum, and StableRank** … we observed a **positive correlation
  > between unsupervised metrics and supervised performance, supporting the utility of unsupervised
  > metrics as effective proxies for model quality assessment**."

Peer-reviewed, genomics, RankMe, concludes it works. This is the single most awkward citation for a
paper claiming rank fails on molecular representations. **Its exploitable weakness: it compares six
*different* pre-trained models — cross-method, which is exactly the regime RankMe itself disclaims
("should only be used to compare different runs of a given method"). Our comparison is within-method.
Say so explicitly.**

### C3. A matched-arm design where effective rank tracks transfer

**Ruan, Zhang, Wang & Zhang, "Muon Learns More Robust and Transferable Features than Adam"**,
arXiv:2606.09658 (S2 CorpusId 289099544). `[NOT INDEPENDENTLY VERIFIED BY ME — S2 record only.]`

> "This transferability advantage is further supported by the diversity of hidden states across
> layers, **as measured by effective rank** … we **prove** that Muon attains larger margins and
> **higher effective rank** than Adam and SGD."

Same architecture, optimiser varied — a matched-arm design in which rank tracks quality, with a proof.

### C4–C7, further defences worth citing (all `[S2-record only, not independently verified]`)

- **Zhang, Deidda, Higham & Tudisco**, "Are We Measuring Oversmoothing in Graph Neural Networks
  Correctly?", arXiv:2502.04591 — *"**rank-based metrics consistently capture oversmoothing, whereas
  energy-based metrics often fail** … drops in the rank align closely with performance degradation"*.
  Rank beating its alternatives, with theory.
- **Zhuo, Wang, Ma & Wang**, "Towards a Unified Theoretical Understanding of Non-contrastive Learning
  via Rank Differential Mechanism", ICLR 2023, arXiv:2303.02387 — the strongest *theoretical* case
  that rank is causally tied to representation quality: *"This rank difference will **provably lead to
  an improvement of effective dimensionality**"*.
- **Kim, Kokilepersaud, Prabhushankar & AlRegib**, "Countering Multi-modal Representation Collapse
  through Rank-targeted Fusion", **WACV 2026**, DOI `10.1109/WACV61042.2026.00461`, arXiv:2511.06450 —
  peer-reviewed, rank used prescriptively, +3.74% SOTA. Relevant because ours is also a multi-modal
  fusion setting.
- **Sun, Lin, Zhang, Duan & Liu**, "Local Dimension Enhancement Representation Learning for
  Skeleton-Based Action Segmentation", **IEEE TIP 2026**, DOI `10.1109/TIP.2026.3682105` — effective
  rank as both diagnostic and training objective, peer-reviewed journal.
- **Billa**, "The Geometric Anatomy of Capability Acquisition in Transformers", arXiv:2602.15997 —
  *"of the geometric measures tested, **only rankme reliably precedes capability acquisition for hard
  tasks**"*, and RankMe is the survivor among competing geometric measures. Contains a limit we can
  use: *"For easy tasks … no precursor is detectable."*
- **Gupta**, "The Geometry of Saturation: Effective Rank Predicts When Labels Stop Helping in Few-Shot
  Classification", arXiv:2606.24903 — effective rank as a *stopping rule that works*, with statistics
  (ρ_pool = 0.6366, p = 2.9×10⁻⁵⁷; AUC 0.787). `[Unrefereed preprint.]`

---

## 4. PARTIAL OVERLAP — must cite (selection; full list in §6)

- **Dai, Xu, Wen, Liu & Huang**, "Exploring Structural Degradation in Dense Representations for
  Self-supervised Learning", arXiv:2510.17299 — within-method **checkpoint-selection** rule built
  from *"a class-relevance measure **and** an effective dimensionality measure"*, +3.0 mIoU over
  end-of-training. Structurally the closest analogue to our contribution outside speech, and its
  design concedes that dimensionality alone is insufficient.
- **Li, Agrawal, Ghosh, Teru, Santoro, Lajoie & Richards**, "Tracing the Representation Geometry of
  Language Models from Pretraining to Post-training", arXiv:2509.23024 — measures **both RankMe and
  α-ReQ** and finds a *"consistent **non-monotonic** sequence of three geometric phases"*, with a
  *"compression-seeking"* phase that **reduces** dimensionality while *"marked with significant
  improvement in downstream task performance"*. **Pre-empts the non-monotonicity limb within a single
  pretraining run**, and note the author overlap with α-ReQ (Agrawal, Ghosh, Richards).
- **Kokilepersaud, Prabhushankar & Al-Regib**, "AdaDim", arXiv:2505.12576 — *"**the best performing
  SSL models do not have the highest H(R) nor the lowest I(R;Z)**, but effectively arrive at a
  balance"*.
- **Jamali, Cheng & Vargas-Hernández**, "Spectral Analysis of Molecular Features: When Richer
  Features Do Not Guarantee Better Generalization", arXiv:2510.14217 — *"**richer spectral features
  do not consistently yield better generalization performance, contradicting common representation
  heuristics used in self-supervised learning**"*, with *"consistently negative correlations"* for
  local 3D representations. **The closest published thing to our claim in a molecular-science
  domain**; cites RankMe and α-ReQ.
- **Yusupov et al.**, "Geometric Metrics and LLMs: What They Measure and When They Work",
  arXiv:2509.25359 — *"apparent discriminative power **collapses once length is controlled**"*. Our
  confound/matched-arm argument, executed for LLMs; strong support for the framing.
- **Adilova, Petzka, Fischer & Geiger**, "Geometric and Information Compression of Representations in
  Deep Learning", arXiv:2606.21593 — *"**low MI does not reliably correspond to geometric
  compression** … a negative and nonlinear relationship **that can reverse when varying training
  setup**"*. A sign reversal under a controlled training change — the same shape of claim as ours, one
  level up.
- **Zaiem, Kemiche, Parcollet, Essid & Ravanelli**, two papers — Interspeech 2023 (arXiv:2306.00452,
  DBLP `conf/interspeech/ZaiemKPER23`) and *Computer Speech and Language* (arXiv:2308.14456):
  *"**altering the downstream architecture structure leads to significant fluctuations in the
  performance ranking of the evaluated models**"*. **Cite defensively** — these attack the *ground
  truth* any rank-vs-downstream correlation is measured against, including ours.
- New-metric suspects, each motivated by an existing proxy being inadequate, each concurrent work:
  **IdEst** (Mordacq, Kalogeiton & Oudot, arXiv:2606.03338 — cites RankMe, LiDAR *and* α-ReQ, the
  densest overlap in the corpus), **Persistence** (Shestov et al., arXiv:2512.15285),
  **Tsitsulin, Munkhoeva & Perozzi** (TAG-ML 2023, arXiv:2305.16562, *"while there is no free
  lunch…"*), **CLID** (Lu et al., TMLR, DBLP `journals/tmlr/LuLBLCS23`), **Q-Score** (Kalibhat et al.,
  AAAI 2024, DOI `10.1609/aaai.v38i12.29201`), **Hart & Tavolara** (arXiv:2407.21590),
  **Darrin et al.** (arXiv:2406.07640, validated in *"natural language processing and molecular
  biology"*).

---

## 5. DOMAIN-ADJACENT — the gap that is genuinely open

Computational pathology and single-cell/omics representation learning are live (Jaume et al., CVPR
2024 `10.1109/CVPR52733.2024.00920` and ECCV 2024 arXiv:2408.02859; HEST-1k NeurIPS 2024
arXiv:2406.16192; DenAdel et al., **Nature Methods 2026**, DOI `10.1038/s41592-026-03120-y` — *"single
-cell foundation models show **no clear data scaling laws**"*; kaiko.ai arXiv:2404.15217; Gustafsson
et al. arXiv:2410.00945 on WSI→gene-expression, closest to our setup).

**Across all 453 de-duplicated citers, none evaluates rank as a selection rule on matched arms in
pathology or transcriptomics.** That gap is real and is where contribution (1) should be positioned.

**One item to check before submission:** Luisto et al., "Domain Fine-Tuning FinBERT on Finnish
Histopathological Reports: Train-Time Signals and Downstream Correlations", arXiv:2604.14815 —
*"predict the benefit of domain-specific pre-training … from observing the **geometry of embedding
changes**"*. Pathology **plus** a rank-geometry selection question. Small and unrefereed, but it is
the one item in the sweep that could scoop the framing and it was triaged at abstract level only.
**Read the full text.**

---

## 6. Required changes to draft §2.2 (for whoever owns the draft)

1. **Delete or rewrite** the sentence "Every one of these tests rank *across* methods, *across*
   checkpoints, or in domains other than joint-embedding representation selection. None tests the
   within-method, matched-arm, fixed-architecture regime that RankMe reserves for itself". It is
   falsified by A1 and strained by A2 and by LiDAR's own hyperparameter claim.
2. **Promote A1 (Aldeneh et al., ICASSP 2025) to the head of §2.2** and state plainly what it
   pre-empts, including the low-rank/high-information limb.
3. **Add a defence subsection** carrying C1 (Deng et al.), C2 (Awasthi et al., BMC Genomics) and C3
   (Ruan et al.), and rebut each. C2 is in-domain and peer-reviewed; it cannot be left out.
4. **Re-scope the novelty claim** to: within-method matched-arm evaluation of rank as a selection
   rule **in computational pathology / transcriptomics**, with a measured seed-level reproducibility
   floor on the statistic and a variance decomposition separating the arm term from the seed term.
   Both of those remain unmatched in the corpus.
5. The `[SEARCH INCOMPLETE]` flag can be **downgraded, not removed**: the citation-graph leg is now
   closed for RankMe, Roy & Vetterli, LiDAR and α-ReQ, but the sweep is abstract-level (12 of 159
   RankMe citers had no abstract at all) and OpenAlex remains an unusable second source. Record the
   residual risk rather than declaring the search complete.

## 7. Provenance

Endpoints used: `https://api.semanticscholar.org/graph/v1/paper/CorpusID:<id>/citations?fields=title,abstract,year,venue,externalIds,authors&limit=1000`,
`https://api.openalex.org/works?filter=cites:<W...>`, `https://api.crossref.org/works/<doi>`,
`https://export.arxiv.org/api/query?id_list=<id>`, `https://dblp.org/search/publ/api`.

Items **I verified personally, from source, rather than relying on the sweep**: A1 (Crossref + arXiv
Atom), A3 (arXiv Atom), C1 (arXiv Atom), C2 (Crossref). Every other item is marked with the record it
rests on. Nothing in this entry is recalled from memory.
