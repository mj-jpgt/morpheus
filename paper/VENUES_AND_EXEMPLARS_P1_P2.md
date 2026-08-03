# Submission venues and structural exemplars — P1 (CALIBRA) and P2 (effective rank)

Compiled 2026-08-03. Scope: venue selection and *structural* exemplar papers to model section order,
framing and phrasing on.

---

## 0. Verification protocol for this document

Three fabricated citations have previously contaminated this project. This document was therefore
compiled under the following rule: **no citation appears here unless it was retrieved from a live
bibliographic service or publisher during compilation.**

Every paper below carries a `VERIFIED` block naming the service used and listing exactly which
fields were confirmed. Where a field could not be retrieved it is marked `NOT VERIFIED` rather than
filled in from memory. Two guesses made during compilation were caught and corrected by this process
and are recorded in §7 as evidence that the protocol is live, not decorative.

Services used, all public and unauthenticated:

* Crossref REST API — `https://api.crossref.org/works/{DOI}` and `?query.title=`
* Europe PMC REST API — `https://www.ebi.ac.uk/europepmc/webservices/rest/search` and `/{PMCID}/fullTextXML`
* NCBI E-utilities / PMC — `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc`
* arXiv API — `https://export.arxiv.org/api/query`
* DBLP — `https://dblp.org/search/publ/api`
* OpenReview API — `https://api.openreview.net/notes/search`
* ar5iv rendered full text — `https://ar5iv.labs.arxiv.org/html/{arXiv-id}`
* PMC article pages — `https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/`

Section orders, figure counts and word counts below were computed by parsing the retrieved full text
(JATS XML section `<title>` elements, or rendered `<h1>`–`<h6>`), not by reading a PDF by eye. Word
counts are **approximate** and are stated with their basis (e.g. "body text including Methods and
figure captions"). Quoted phrases are verbatim from the retrieved full text.

**What is NOT verified anywhere in this document:** any statement about how a venue will *receive*
these papers. Editorial-policy quotes and published precedents are evidence; they are not
predictions.

---

## 1. Lead with the bad news — placement risk

Read this before the ranked lists.

### 1.1 P1 has a structural placement problem that is not about quality

P1 is a *corrective* paper with **no external cohort**. In computational pathology the external
validation cohort has become a de facto admission requirement, and — importantly — the exemplar
literature says so out loud. Howard et al. 2021, the closest published analogue to P1, devotes a
passage of its Discussion to exactly this bind:

> "When developing predictive histologic models for a large number of features, external validation
> of every finding can be impractical/infeasible. Furthermore, adequate external validation datasets
> may not be readily available for rare cancer subtypes."

(Verbatim, Howard et al. 2021 Discussion, retrieved from Europe PMC full text — see §3.1.)

That is the single most useful sentence found during this search. It is a published, peer-reviewed,
high-profile precedent for the exact limitation P1 must declare, written by authors in P1's own
subfield, and it survived review at Nature Communications. **P1 should cite it at the point where it
declares the scope limitation.** It converts "we lack external validation" from an apology into an
argued position with a citation behind it.

But the risk is not only a reviewer risk. **At one venue it is a written rule.** *Bioinformatics*
(OUP) states: *"Machine learning papers **must** report the performance on an independent test set.
It is not sufficient to report the average error over the individual cross-validation sets."*
(<https://academic.oup.com/bioinformatics/pages/author-guidelines>). P1 would be in violation on
submission. Two further gate types disqualify P1 elsewhere — a "clear advance over state of the art"
mandate (*Medical Image Analysis*, *Genome Biology* Methodology) and a clinical-validation scope
(*npj Digital Medicine*). All three are quoted and sourced in §2.0.

The good news buried in that: **no venue examined rejects negative results as such.** The searches
found no anti-negativity policy anywhere. What they found were novelty and validation gates, which is
a different problem with a different fix — pick venues without those gates (§2), and argue the
*instrument* is novel even though the *result* is deflationary.

### 1.2 P1 has an unlisted prior-art exposure that should be checked before anything else

The venue search turned up a paper that is not in `P1_CALIBRA_DRAFT.md`'s reference list and is the
closest published relative of CALIBRA found anywhere:

**Hamdan, Love, von Polier, Weis, Schwender, Eickhoff & Patil, "Confound-leakage: confound removal in
machine learning leads to leakage", *GigaScience* 12: giad071.**
DOI <https://doi.org/10.1093/gigascience/giad071> (verified via Europe PMC, PMC10541796, CC BY, full
text retrieved).

Its abstract states that "contrary to general expectations, linear confound regression can **increase**
the risk of confounding when combined with nonlinear ML approaches… information leaked via CR can
increase null or moderate effects to near-perfect prediction." That is adjacent to, and possibly
overlapping with, draft §2.4 (residualisation-induced correlation), which is one of P1's stronger
novelty claims. **Read it before finalising §2.4.** Details and its retrieved section order are in
§2.2.

### 1.3 P2 has a systematic venue-class problem

Three concrete findings, in descending order of how much they should change your plans.

**(a) One venue has a written rule that disqualifies P2 outright, and it is the one you might have
assumed was safest.** JMLR states: *"JMLR publishes papers on the theory and methods of machine
learning but **does not publish applications of machine learning to other domains**. JMLR favors
papers of interest to a broader machine learning audience and may deem a paper unsuitable if the
editorial board finds its audience too narrow."*
(<https://www.jmlr.org/author-info.html>) P2's evidence base is a histology/transcriptomics
application. The same page notes the editor-in-chief can reject "without written review". Do not
submit there. Note the shape of this risk: it is not that the result is negative, it is that the
*demonstration is applied*. That risk follows P2 to any venue that prizes generality, and it is the
one to write defensively against.

**(b) Conferences did not become safe for negative results — they became *explicit* about them, and
raised the bar.** NeurIPS added a formal "Negative Results" contribution type for 2026, which is
genuine progress. But the reviewing rubric says, twice: "The significance bar for Negative Results
papers is high" and "The originality bar for Negative Results papers is high" — an elevated bar
stated for only two of five contribution types.
(<https://neurips.cc/Conferences/2026/ReviewerGuidelines>) Meanwhile ICML requires reviewers to score
**Originality and Significance separately on 1–4 scales**
(<https://icml.cc/Conferences/2026/ReviewerInstructions>), which is the concrete mechanism by which a
correct-but-unflattering result scores "2: fair" and dies. The documentary contrast with TMLR, which
forbids novelty from being grounds for rejection at all, is the whole basis of the ranking in §4.

**(c) There is a closer published precedent than expected — and it is a competitor, not just a
citation.** Fang, Li, Sun & Wang, "Rethinking the Uniformity Metric in Self-Supervised Learning",
ICLR 2024 (<https://openreview.net/forum?id=3pf2hEdu8B>, arXiv:2403.00642; *verified independently
via DBLP: venue ICLR, year 2024, four authors*). It shows that the widely used uniformity metric fails
to satisfy desirable properties and mishandles dimensional collapse. This both proves the genre is
placeable at a top conference and means P2 must differentiate itself explicitly — see §5.0.

The known precedents (Leavitt & Morcos at ICLR 2021; Kynkäänniemi et al. at ICLR 2023; Musgrave et al.
at ECCV 2020) are real evidence that the genre gets in. They are, however, largely from authors at
large industrial labs, which is a confound in the precedent itself and is stated here rather than
glossed.

### 1.4 P1 and P2 currently collide, and it is a dual-submission problem

This was not in the brief but it is the most actionable thing found, so it goes near the top.

`paper/P1_CALIBRA_DRAFT.md` §4.11 is titled "Effective rank does not track information content — four
independent instances" and contains a four-row table of exactly the dissociations that constitute P2
(rank +107% / specificity flat; rank −17% / channel flat; rank pinned at 16/16 while retrieval falls
to 0.000 below a chance of 0.062; rank −18% / channel −67%), plus the seed-instability figures
(9.14 → 34.12). P2 is that table, expanded.

P1 §4.11 already carries a disclaimer — "A fuller treatment belongs to a companion paper and is not
claimed here" — which handles the *claim* but not the *reuse*. TMLR's rule is about results, not
claims. Verbatim from <https://jmlr.org/tmlr/editorial-policies.html>:

> "Unlike many other journals, TMLR only accepts original contributions that don't reuse the authors'
> own prior work. In particular, we do not accept submissions that are expanded versions of conference
> papers. There should not be any reuse of written text, figures or **results** between the submitted
> paper and any paper which has been published, accepted for publication, **or submitted in parallel
> at another archival, peer-reviewed venue**."

So: if P1 is under review at any archival peer-reviewed journal while P2 is under review at TMLR, and
both contain the four-instance table, that is a dual-submission violation on its face. The same clause
also confirms the safe harbour — "It is acceptable for a submission to overlap with the author's
previous work if it was shared at venues or tracks that are publicly declared, in writing, to be
non-archival, such as workshops, or on preprint servers such as arXiv and bioRxiv" — so **arXiv
preprints of either paper are not a problem.**

Three ways out, in order of preference:

1. **Cut §4.11 from P1 down to a single forward-reference sentence** with no table and no numbers,
   citing the P2 preprint. P1 loses nothing it needs: §4.11's own text concedes the point is "narrow"
   for P1's purposes.
2. **Sequence them** — P2 to TMLR first, P1 submitted after P2 is public, with §4.11 reduced to a
   citation.
3. **Send P2 somewhere without that clause** (ICLR, NeurIPS E&D). This is the worst option because it
   gives up the single strongest reason to prefer TMLR (§4).

Resolve this before either paper goes out. It is cheap now and expensive later.

### 1.5 One finding that is good news for P2

A targeted search of arXiv (abstract-level, via the arXiv API) for an existing paper that debunks
effective rank as a representation-quality proxy returned **nothing**. Queries run:
`all:"effective rank" AND all:"self-supervised"`, `ti:"effective rank" AND abs:"representation"`,
`abs:"RankMe"`, `abs:"effective rank" AND abs:"does not"`,
`abs:"representation quality" AND abs:"proxy" AND abs:"rank"`. The returned hits were all papers that
*use* effective rank as a diagnostic or *regularise* toward it — none that dissociate it from
information content.

Caveat, stated plainly: arXiv abstract search is not exhaustive, it does not cover non-arXiv venues,
and the search was constrained (the WebSearch budget for this session was exhausted, so no
general-purpose web search was available). Treat this as "no prior art found by these five queries",
not "no prior art exists". A proper novelty sweep with full-text search should be run before
submission.

---

## 2. P1 — venues

**Venue class verdict: computational-biology / data-science journals with a *resource* or
*re-evaluation* article type — not the computational-pathology imaging mainstream, and not the
bioinformatics mainstream either.** The reason is specific and evidenced below: no venue examined
rejects negative results as such, but three separate gate types will kill P1 — an explicit
**independent-test-set mandate**, a **"clear advance over state of the art" mandate**, and a
**clinical-validation scope**. P1 fails all three by construction, and they are unevenly distributed
across venues.

> **Verification note.** Policy quotations below were retrieved from the URLs given. Where a figure
> came from a third party (DOAJ) or an archived snapshot rather than the publisher's live page, that
> is stated inline. Publisher sites that block automated retrieval (cell.com, sciencedirect.com,
> modernpathology.org) are flagged where they affect confidence.

### 2.0 Read this first — the three gates that disqualify P1, with the venues that impose them

**Gate 1 — independent test set required. This is the single hardest disqualifier found.**
*Bioinformatics* (OUP) states, verbatim:

> "Machine learning papers must report the performance on an independent test set. It is not
> sufficient to report the average error over the individual cross-validation sets."
> <https://academic.oup.com/bioinformatics/pages/author-guidelines>

P1 has no external cohort and says so. Do not submit to *Bioinformatics*. This is the one venue where
the declared scope limitation is a written rule violation rather than a reviewer risk.

**Gate 2 — must demonstrate a clear advance over prior methods.** A corrective instrument has nothing
to beat, so this gate is structurally hostile:

* *Medical Image Analysis*: "Before submitting a manuscript to Medical Image Analysis, please ensure
  that your submission contains a 'significant methodological contribution'… **Manuscripts which do
  not contain a significant methodological contribution will be rejected without review.**"
  (Guide for Authors; retrieved from an Internet Archive snapshot dated 2024-02-24 —
  <http://web.archive.org/web/20240224213732/https://www.sciencedirect.com/journal/medical-image-analysis/publish/guide-for-authors>
  — because sciencedirect.com blocks automated retrieval. **Currency of this quote is NOT VERIFIED
  against the live page.**)
* *Genome Biology*, Methodology type: "Methodology articles should describe novel methods that are
  shown to be a **clear advance over existing state-of-the-art methods in a side-by-side
  demonstration**… For computational methods, demonstration of superiority should be carried out
  using the same dataset."
  <https://link.springer.com/journal/13059/submission-guidelines/methodology>
* *PLOS Computational Biology*, Methods/Software types: "Enhancements to existing published methods
  or software will only be considered if those enhancements bring exceptional new capabilities."
  <https://journals.plos.org/ploscompbiol/s/journal-information>

**Gate 3 — clinical application and validation.** *npj Digital Medicine*: "The journal typically does
not consider clinical research using off-the-shelf digital tools and artificial intelligence models,
pre-clinical basic science studies, purely observational studies, case studies, and small-scale
preliminary studies", with scope bullets framed as "Clinical application of novel and **validated**
artificial intelligence and machine learning models."
<https://www.nature.com/npjdigitalmed/aims> — despite the journal having published exactly P1's genre
in 2019 (Badgeley et al., below). The 2019 precedent does not survive today's scope statement.

**A structural trap worth knowing now.** "Matters Arising" and equivalent comment formats are useless
to P1, because they are restricted to critiquing the *host journal's own* papers: "Nature
Communications does not consider Matters Arising on papers published in other journals"
(<https://www.nature.com/ncomms/submit/matters-arising>); "Genome Biology does not consider Matters
Arising on papers published in other journals"
(<https://link.springer.com/journal/13059/submission-guidelines/matters-arising>). P1 critiques a
practice spread across many journals, so it must enter as a **primary research / resource article**
everywhere.

**And the headline: no venue examined has a policy rejecting negative results as such.** Only two
venues have affirmative written policy in P1's favour — GigaScience's anti-impact criterion and
eLife's Replication Studies clause. Everywhere else, the evidence is published precedent.

### 2.1 Rank 1 — **Nature Communications**

**Scope fit (judgement).** Best combination of domain match, format capacity and precedent. Scope is
"all areas of the biological, health, physical, chemical, Earth, social, mathematical, applied, and
engineering sciences" (<https://www.nature.com/ncomms/aims>).

**Negative/corrective work — two verified in-domain precedents.**

* **Howard et al. 2021** (§3.1) — the site-confound paper. A purely corrective paper about TCGA-based
  deep-learning pathology, with P1's exact limitation, published here.
* **Chuhan Wang, Adam S. Chan, Xiaohang Fu, Shila Ghazanfar, Jinman Kim, Ellis Patrick et al.,
  "Benchmarking the translational potential of spatial gene expression prediction from histology",
  *Nature Communications* 16, article 1544 (2025).**
  DOI <https://doi.org/10.1038/s41467-025-56618-y>
  *VERIFIED by me via Crossref: title, journal, volume 16, article number 1544, date 2025-02-11.*
  This benchmarks eleven methods on P1's exact task family and should be cited in P1's Related Work
  regardless of venue.

Policy support exists but is narrow: "The journal is particularly interested in publishing clinical
trials, including trials with **negative results** of high relevance to a particular community"
(<https://www.nature.com/ncomms/submit/clinical-research>) — that is the clinical-trials section, not
a general negative-results policy. Do not over-read it.

**Article type and limits** (<https://www.nature.com/ncomms/submit/article>), verbatim: main text
"should be limited to **5,000 words**" excluding Abstract, Methods, References and Figure legends;
abstract "no more than 200 words"; "Methods are typically less than **3000 words**"; "As a guide,
references should not exceed **70**"; "Articles may have up to **10 display items**". The journal "is
flexible with regard to the format of initial submissions" and "does not consider pre-submission
enquiries."

**This envelope is the reason it ranks first:** 5,000 words *plus* a 3,000-word Methods *plus* ten
display items is more capacity than any other high-prestige venue on this list, and it is what Howard
et al. used (≈10,900 body words including Methods and legends, 6 figures, 80 refs).

**OA / APC.** Fully open access. "The current APC, subject to VAT or local taxes where applicable, is:
**£5490.00/$7350.00/€6150.00**" — <https://www.nature.com/ncomms/about/article-processing-charges>
(retrieved 2026-08-03). CC BY or CC BY-NC-ND.

**Time to first decision.** **9 days median** submission → first editorial decision; 242 days median
submission → acceptance (<https://www.nature.com/ncomms/journal-impact>). *Caveat: the journal's own
/about page states 8 days; the discrepancy is on the publisher's site and is unresolved.* Howard et
al.'s own record (received 2020-12-09, accepted 2021-07-01 ≈ 6.7 months) is consistent with the 242-day
figure.

**Flags.** A soft novelty gate — "Papers published by the journal aim to represent important advances
of significance to specialists within each field" (<https://www.nature.com/ncomms/aims>). **No
external-validation requirement and no clinical-utility requirement were found.** The real constraint
is that the full must-PASS/must-FAIL battery will have to live in Supplementary.

### 2.2 Rank 2 — **GigaScience** — and it has published P1's nearest relative

**This is the most important discovery in the P1 venue search, on two counts.**

**Count one — it is the only venue with an affirmative written policy removing the impact gate.**
Verbatim:

> "Criteria for publication are **reproducibility, usability and utility, rather than subjective
> assessment of 'impact'**." — <https://academic.oup.com/gigascience/pages/About>

That is the closest thing in the journal world to TMLR's protection for P2 (§4).

**Count two — it published the closest published relative of CALIBRA anywhere, and P1 must cite it.**

**Sami Hamdan, Bradley C. Love, Georg G. von Polier, Susanne Weis, Holger Schwender, Simon B.
Eickhoff, Kaustubh R. Patil. "Confound-leakage: confound removal in machine learning leads to
leakage." *GigaScience* **12**: giad071.**
DOI <https://doi.org/10.1093/gigascience/giad071> · PMID 37776368 · PMCID PMC10541796

> **VERIFIED by me** via Europe PMC (PMC10541796, `isOpenAccess: Y`, CC BY; full text retrieved).
> Confirmed: title, all seven authors, journal, volume 12, DOI, PMID, PMCID.
> *Year ambiguity:* the Europe PMC `pubYear` field reads **2022** while the article ID (`giad071`) and
> PMID indexing indicate **2023**. Resolve before citing — do not guess.
> Retrieved section order: Background / Results / Conclusions (structured abstract) → Introduction →
> Results (Walk-through analysis; TaCo removal for binary classification; Confound removal for
> regression; Analyses of benchmark data; TaCo removal increases performance of nonlinear methods;
> CR using weaker confounds also increases performance; Increased performance after TaCo removal is
> due to confound-leakage; Possible mechanisms for confound-leakage; …deviation from normal
> distributions; …limited precision features; Confound-leakage poses danger in clinical
> applications) → Discussion → Conclusions and future directions → Methods.

**Why this matters beyond venue choice.** Its finding, from the retrieved abstract, is that "contrary
to general expectations, linear confound regression can **increase** the risk of confounding when
combined with nonlinear ML approaches… information leaked via CR can increase null or moderate effects
to near-perfect prediction." That is adjacent to — and possibly overlapping with — P1's §2.4 on
residualisation-induced correlation. **This is a potential prior-art exposure and it is currently
absent from the draft's reference list.** Read it before finalising §2.4's novelty claim. Its Results
ordering (walk-through → mechanism → benchmark replication → mechanisms → danger) is also a usable
structural template.

A second relevant GigaScience paper: Tamas Spisak, "Statistical quantification of confounding bias in
machine learning models", *GigaScience* **11** (2022), DOI
<https://doi.org/10.1093/gigascience/giac082>. *Citation reported by the venue-research agent from
Crossref; not independently re-verified by me.*

**Article types and limits** (<https://academic.oup.com/gigascience/pages/Instructions_To_Authors>).
**Research**: abstract 250 words max, **no main-text word limit stated**
(<https://academic.oup.com/gigascience/pages/research>). **Technical Note**: same, and the natural
home for the released instrument. Commentary: max 1,200 words, max 10 references.

**OA / APC.** Fully OA, CC BY 4.0. **≈$2,638 USD — DOAJ record only** (ISSN 2047-217X, record updated
2026-01-15). The OUP charges page renders prices through a JavaScript widget and returned a
placeholder, so the **publisher-confirmed figure is NOT VERIFIED**. This is nonetheless the cheapest
option on the list by a wide margin.

**Time to first decision. NOT VERIFIED** — no publisher metric exists. Received→accepted proxies from
article records gave medians of ≈129 days (n=4) and ≈202 days (n=9); both measure the full cycle, not
first decision.

**Flag — one real trap.** "Data validation, basic comparisons to other data, etc. are suitable for a
Data Note, but are not alone sufficient for a Research Article", and Research Articles must "provide
some scientific insight and conclusions"
(<https://academic.oup.com/gigascience/pages/research>). **Mitigation:** frame the 0.147 chance
baseline and the 76–82% random-gene-set result as *findings about the field*, not as validation of an
instrument. No "significant advance" gate, no novelty gate and no external-validation requirement were
found anywhere in GigaScience's policies.

### 2.3 Rank 3 — **Cell Reports Methods**, submitted as a **Resource**

**The article-type definition is close to bespoke for P1.** Verbatim:

> "**Resources** at Cell Reports Methods are empirical studies that **inform the use and
> interpretations of existing methods**. They **do not present a new method per se**, but provide
> authors with a forum in which they may **re-evaluate current approaches in a data-driven way**."
> <https://www.cell.com/cell-reports-methods/information-for-authors/article-types>

**Precedent.** **David Wissel, Daniel Rowson, Valentina Boeva, "Systematic comparison of multi-omics
survival models reveals a widespread lack of noise resistance", *Cell Reports Methods* **3**(4):100461
(2023), DOI <https://doi.org/10.1016/j.crmeth.2023.100461>.** *VERIFIED by me via Crossref: title,
three authors, journal, volume 3, article 100461, 2023.* TCGA-based, and asks P1's question — does the
extra signal actually add anything.

**Limits.** Article and Resource: "under 7,000 words, not including the reference list or the STAR
Methods", "no more than 7 figures and/or tables". **Mandatory limitations section:** "We require a
'limitations of the study' subsection at the end of the discussion" — P1's no-external-cohort
limitation has a designated, expected home here, which is a real advantage.

**APC.** **$4,970**, "excluding taxes" (<https://www.cell.com/open-access>). *Currency inferred as USD
from Elsevier's general pricing statement, not from a per-journal statement — see §7.*

**Speed.** 7 days median to first editorial decision (desk). Post-review median not published.

**Flag.** The **Article** type carries a validation clause ("the presentation of the method, a
comparison of the method with related methods when appropriate, and validation"). The **Resource**
description does not. Submit as a Resource.

### 2.4 Rank 4 — **Patterns** (Cell Press)

**Precedent — the field's canonical example of this genre.** Sayash Kapoor & Arvind Narayanan,
"Leakage and the reproducibility crisis in machine-learning-based science", *Patterns* **4**(9):100804
(2023), DOI <https://doi.org/10.1016/j.patter.2023.100804> (PMID 37720327, PMCID PMC10499856).
Published as a full research-format article; finds leakage across 17 fields and 294 papers.

**Limits.** Research article "approximately 5,000–7,000 words" with "**There is no strict limit on the
number of figures, tables, or references**"
(<https://www.cell.com/patterns/information-for-authors/article-types>). Mandatory "bigger picture"
box, max 300 words. Replication policy: "Although Cell Press journals do not have a dedicated format
for replication studies, we are open to considering them."

**APC.** **$4,900**, "excluding taxes" (<https://www.cell.com/open-access>).

**Speed.** 7 days median to first editorial decision; the FAQ states "The first decision after peer
review is usually rendered within 45 days."

**Flag — selectivity is the risk.** "The journal is highly selective and seeks to publish works that
offer major advances", and suitability is judged on "the conceptual novelty of the work". Validation
language is soft ("The methods presented in our papers are expected to be well validated") with **no
external-cohort requirement found**.

### 2.5 Rank 5 — **PLOS Computational Biology**

**Two decisive facts, one good and one bad, and both are commonly got wrong.**

**Good — it neutralises P1's biggest weakness in writing:**

> "Inclusion of **experimental validation is not required for publication**, but should be referenced
> where possible." — <https://journals.plos.org/ploscompbiol/s/journal-information>

**Bad — the famous PLOS negative-results clause is *PLOS ONE's*, not this journal's.** PLOS ONE:
"Studies reporting negative results — In keeping with our mission to publish all valid research, we
consider negative and null results"
(<https://journals.plos.org/plosone/s/criteria-for-publication>). *PLOS Computational Biology* has no
such clause; its stated Criteria for Publication are "**Originality; Innovation; High importance to
researchers in the field; Significant biological and/or methodological insight; Rigorous methodology;
Substantial evidence for its conclusions**" (journal-information URL above, verified independently by
me). Do not assume the PLOS ONE language transfers — it does not.

**Format — the most permissive on this list.** Research Article: "**Manuscripts can be any length.
There are no restrictions on word count, number of figures, or amount of supporting information.**"
Abstract ≤300 words. Note the Benchmarking type is *not* usable: "Novel tools created by the authors
in parallel with the benchmark should not be included in the article."

**Precedent.** Venet et al. 2011 (§3.2) is here, which is direct evidence the random-control argument
publishes at this venue at full length (received 2011-04-27, accepted 2011-09-07 ≈ 4.4 months). A
recent PCB *Research Article* whose central contribution is a confound correction: **NOT VERIFIED** —
the recent corrective item found (Chicco & Jurman 2025, DOI 10.1371/journal.pcbi.1013673) is a
**Formal Comment**, not a Research Article.

**APC.** **$3,165 USD** for a PCB Research Article (<https://plos.org/publish/fees/>).

**Speed.** "Time to first decision **41 days**", time to publication 230 days, acceptance rate 33%
(2025 figures, <https://plos.org/metrics/>). This is the fastest *verified, publisher-stated*
first-decision figure that is not a desk-screen number.

**Flag.** "Innovation" and "High importance" are explicit gates, and "biological significance" is a
named criterion — a pure instrument paper with no biological finding is exposed on all three.

### 2.6 Rank 6 — **Nature Methods**

**Superb article-type fit, poor aims-and-scope fit — a high-variance submission.**

**In favour.** From <https://www.nature.com/nmeth/content> (retrieved independently by me; the page is
static):

* **Brief Communication** — "a concise report describing potentially groundbreaking yet preliminary
  method or tool developments, highly practical tweaks to an existing method or tool, software
  platforms, resources of broad interest, and **technical critiques of widely used methodologies**."
  Abstract ≤70 w; main text 1,200 w (up to 1,600 at discretion) **including** abstract, refs and
  legends; max 2 display items (3 at discretion); ~20 refs.
* **Analysis** — "a report presenting comprehensive performance comparisons of established, related
  methods or tools, of key importance to a field of research." Abstract ≤150 w; 3,000 w (up to 5,000
  at discretion) *excluding* abstract, Methods, refs and legends; ≤6 display items; ~50 refs.
* **Article** — "a report describing a novel method or tool… include strong validation data to
  demonstrate performance, reproducibility, general applicability". Same envelope as Analysis.
* **Registered Report** — "Stage 2 manuscripts containing Results and Discussion are also peer-reviewed
  but **cannot be rejected for reasons of novelty or perceived importance of the results**." *Not
  available to P1* — the same page states completed comparative analyses "should be submitted as
  Analyses" — but it is the right vehicle for any future CALIBRA-based study not yet run, and it is
  the journal-world analogue of TMLR's protection.

Up to 10 Extended Data items are permitted
(<https://www.nature.com/nmeth/submission-guidelines/preparing-your-submission>), which materially
relieves the display-item cap.

**Precedent.** Constantin Ahlmann-Eltze, Wolfgang Huber, Simon Anders, "Deep-learning-based gene
perturbation effect prediction does not yet outperform simple linear baselines", *Nature Methods*
**22**, 1657–1661 (2025), DOI <https://doi.org/10.1038/s41592-025-02772-6> — a Brief Communication
with "does not" in the title. *Verified via Crossref by the venue agent and independently corroborated
in the sibling repo analysis.*

**Flag — the biggest single scope problem in the set.** Aims and scope: "To enhance the practical
relevance of each paper, **the description of the method must be accompanied by strong validation, an
application to an important biological question and results illustrating its performance in comparison
to available approaches**" (<https://www.nature.com/nmeth/aims>). That sentence is written for
method-development Articles, and the Brief Communication "technical critique" route arguably escapes
it — but P1 would be relying on an editor reading it that way.

**APC.** Hybrid. Gold OA **£9,390 / $12,850 / €10,850**
(<https://www.nature.com/nmeth/submission-guidelines/publishing-options>); the subscription route is
free to authors. **Speed:** 10 days median to first decision, 276 days to acceptance
(<https://www.nature.com/nmeth/journal-impact>).

### 2.7 Rank 7 — **eLife** — the structural hedge

**Strongest written tolerance for corrective work, and a guaranteed outcome once reviewed.** Verbatim
from <https://elifesciences.org/about/peer-review>:

> "Publishing with eLife is different because **there is no accept or reject decision after peer
> review**: rather, every article we review is published as a Reviewed Preprint."

> "Once an article has been selected for peer review, the authors can be sure that it will be
> published."

And from the Replication Studies instructions: "Replication is an important part of scientific
progress. We welcome Replication Studies that provide new insights into previously published research
in eLife or other journals… **The outcome of the replication will not affect whether or not we proceed
with peer review.**"

**Best-fit type: Tools and Resources** — "do not have to report major new biological insights or
mechanisms, but it must be clear that they will enable such advances to take place… the new method
must be thoroughly compared and benchmarked against existing methods used in the field." Research
Article: no maximum length, suggested ≤5,000 words excluding Methods/refs/legends, **no limit on
display items**.

**APC.** "$3750 and is payable at this stage (although the fee is **waived for anyone who cannot
afford to pay**)". **Speed: NOT VERIFIED** — no official median; two article-timeline proxies gave 81
and 141 days from sent-for-review to Reviewed Preprint v1.

**Read the gate correctly, and weigh the assessment risk.** The decision has moved to the pre-review
screen — "we will not review papers that appear to report incremental results" — and what publishes
alongside P1 is a permanent **eLife Assessment**. The significance vocabulary (landmark / fundamental
/ important / valuable / useful) is defined entirely in terms of the *implications of findings*, with
no vocabulary for "corrects the field's error", so expect an asymmetric result: a modest significance
word attached permanently, paired with a potentially strong evidence word — "compelling" is defined as
"evidence that features methods, data and analyses more rigorous than the current state-of-the-art",
which a full negative-control battery fits well.
(<https://elifesciences.org/about/elife-assessments>) Note also that "convincing" is pegged to
"appropriate and **validated** methodology", which reviewers may read as licence to demand external
validation.

### 2.8 Wildcard worth a serious look — **PLOS Biology, "Meta-Research Article"**

> "data-driven articles that examine significant questions of **how biological research is designed,
> carried out, communicated and evaluated**"
> <https://journals.plos.org/plosbiology/s/journal-information>

This is the only article type found anywhere whose *definition* is what CALIBRA is. Against it: "PLOS
Biology is highly selective and publishes significant advances resulting from original lines of
inquiry that have a broad impact", and the APC is **$5,500** non-member
(<https://plos.org/publish/fees/>). A long shot, but a legitimate reframing — and the reframing itself
("this is meta-research about measurement practice, not a pathology method") is worth having in the
cover letter wherever P1 goes.

### 2.9 Venues to avoid, and why

| Venue | Disqualifier |
|---|---|
| **Bioinformatics (OUP)** | "Machine learning papers **must** report the performance on an independent test set." Also: papers using leave-one-out are editorially rejected absent special circumstance; ~30% acceptance; 5,000-word cap. |
| **npj Digital Medicine** | Scope is clinical application of *validated* models; explicitly excludes small-scale preliminary and purely observational work. (Its own 2019 precedent — Badgeley et al., *npj Digit Med* **2**, article 31, DOI <https://doi.org/10.1038/s41746-019-0105-1>, *verified by me via Crossref* — is literally P1's chance-baseline argument, but predates the current scope statement.) |
| **Medical Image Analysis** | "Manuscripts which do not contain a significant methodological contribution will be rejected without review." Its published pitfall paper (Rivoir et al. 2024) also delivered a new SOTA — that is the template it rewards. |
| **Genome Biology** | Methodology and Software types require "a clear advance over existing state-of-the-art methods in a side-by-side demonstration"; Research requires "a substantial advance over previous studies". Strong appetite for the genre, wrong gate. |
| **Nature Machine Intelligence** | Article/Analysis capped at 3,500 words / 6 display items / 50 refs; APC £9,390 / $12,850 / €10,850. Excellent precedent (DeGrave et al. 2021), unusable envelope. |
| **Modern Pathology, J. Pathology Informatics** | **NOT ASSESSED** — every guideline page returned 403/Cloudflare to automated retrieval. Named here so their absence is not read as a judgement. |

### 2.10 Recommended sequence for P1

1. **Read Hamdan et al. 2023 (§2.2) before doing anything else.** It is the nearest published relative
   and a possible prior-art exposure for §2.4. Novelty framing depends on what it does and does not
   claim.
2. **Nature Communications first.** Only venue where the format fits without amputation *and* a nearly
   identical corrective paper is in print. Cite Howard et al. at the limitation, using its own
   "external validation of every finding can be impractical/infeasible" sentence.
3. **GigaScience second** — the only affirmative anti-impact criterion, no length cap, cheapest APC,
   and it published the closest relative. Frame the chance baseline and random-gene-set results as
   findings, not validation.
4. **Cell Reports Methods (Resource) third** — the article-type language is close to bespoke and the
   mandatory limitations subsection gives the no-external-cohort declaration a designated home.
5. **Patterns / PLOS Computational Biology** as the next tier; **eLife** as the structural hedge if the
   external-cohort objection proves fatal, accepting that the Assessment travels with the paper.
6. Resolve §1.4 (the P1/P2 result overlap) **before** submitting anywhere.

---

## 3. P1 — structural exemplars

Chosen for **matching kind of contribution**, not topic. P1 has three distinct rhetorical problems
and each exemplar solves one of them better than the others do:

| P1's problem | Exemplar that solves it |
|---|---|
| "A widely reported effect in *this exact subfield* is confounded, and here is the fix" | Howard et al. 2021 |
| "The baseline is wrong, therefore a body of reported numbers is inflated" | Eklund et al. 2016 |
| "Random controls reproduce most of the reported signal" + "the contribution is a negative-control battery" | Venet et al. 2011 |

### 3.1 Howard et al. 2021 — the domain-matched confound exposé with a corrective instrument

**Citation.** Frederick M. Howard, James Dolezal, Sara Kochanny, Jefree Schulte, Heather Chen, Lara
Heij, Dezheng Huo, Rita Nanda, Olufunmilayo I. Olopade, Jakob N. Kather, Nicole Cipriani, Robert L.
Grossman, Alexander T. Pearson. "The impact of site-specific digital histology signatures on deep
learning model accuracy and bias." *Nature Communications* **12**, article 4423 (2021).
DOI: <https://doi.org/10.1038/s41467-021-24698-1>
Full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8292530/>

> **VERIFIED** via Crossref (`/works/10.1038/s41467-021-24698-1`) and Europe PMC (PMC8292530,
> `isOpenAccess: Y`, licence CC BY). Confirmed from source: title, all 13 authors and their order,
> year, journal, volume 12, article number 4423, DOI, PMID 34285218, PMCID.
> Section order, figure count, reference count and word count computed from the retrieved JATS XML.
> Received 2020-12-09; Accepted 2021-07-01 (from the retrieved record).

**Why it is structurally analogous.** This is the closest published relative P1 has. Same repository
(TCGA), same confounder (tissue submitting site), same conclusion shape (models are reading the site,
not the biology), and — critically — the same two-part structure P1 has: *diagnose the confound*,
then *ship an instrument that removes it* (their preserved-site cross-validation via quadratic
programming; P1's cross-fitted residualisation plus certificate). It also deflates named prior claims
without being purely destructive.

**Actual section order** (from JATS `<title>` elements):

1. Abstract
2. Introduction
3. Results
   - Characterization of clinical and digital imaging heterogeneity in TCGA
   - Deep-learning algorithms accurately identify tissue submitting site
   - An artificial simulation of site-specific digital histology signatures
   - Preserved-site cross-validation — a quadratic programming solution
   - Impact of site-specific digital histology signatures on deep-learning model performance
4. Discussion
5. Methods (Patient cohorts; Image processing and deep-learning model; Statistics and reproducibility; Reporting summary)
6. Data availability / Code availability / Competing interests

**Length and figures.** 6 figures, 0 tables, 80 references. Body text ≈10,900 words *including*
Methods and figure captions (parsed from `<body>`). This is a long-format Nature Communications
Article, not a brief report.

**Note the ordering choice, because P1 should copy it.** The corrective instrument
("Preserved-site cross-validation — a quadratic programming solution") appears as the **fourth of
five** Results subsections — after the confound is characterised, after it is shown to be
DL-detectable, and after a synthetic simulation establishes the mechanism. The instrument is
*earned*, not announced. P1's current draft front-loads instrument mechanics; this exemplar argues
for establishing the problem first.

**How it handles our shared rhetorical problem.**

*How it opens:* not with the confound. It opens with clinical grounding — "A standard component of
the diagnosis of nearly all human cancers is the histologic examination of hematoxylin and
eosin-stained tumor biopsy sections" — then walks the reader up through what DL has been claimed to
do, and only then turns. The pivot sentence is:

> "However, the overfitting of digital histology models to site-level characteristics has been
> incompletely characterized and is infrequently accounted for in the internal validation of deep
> learning models."

That is a model for a negative paper's turn: the prior literature is described generously and
accurately for a full paragraph before the problem is named, and the problem is named as an *omission*
("incompletely characterized", "infrequently accounted for") rather than as an error by named
authors.

*How it states the negative claim without overreaching:* it quantifies the scope of the damage and
concedes where the damage is absent. Verbatim:

> "The effect size is small for the majority of features and is absent for most features with a clear
> histologic basis such as tumor histologic subtype and grade."

and

> "Of course, a number of models initially tested in TCGA without preserved-site cross-validation have
> maintained accurate prediction in external validation cohorts, such as prediction of microsatellite
> instability or BRAF mutations in colon cancer."

This is the single most important move for P1 to imitate. Howard et al. explicitly names the cases
where the prior literature *survives* their critique. It costs them nothing and it makes the cases
that don't survive far harder to dismiss as motivated. P1 has an exact analogue available: the
targets whose signal survives residualisation, and the attenuation ≈ 1 result.

*How it presents controls:* the "artificial simulation of site-specific digital histology signatures"
subsection is a planted-signal control — they inject a known synthetic site signature and show the
pipeline recovers it. This is structurally identical to P1's injection certification, and it is placed
*before* the corrective method, as motivation for it.

*How it words its limitations:* the external-validation passage quoted in §1.1, plus a careful hedge
about generalising the mechanism:

> "Given that traditional image and textural characteristics vary between sites in TCGA, it is likely
> that non-deep-learning prognostic studies that predict outcome from traditional image analysis
> features may suffer from a similar bias."

Note "it is likely that … may suffer" — a deliberately double-hedged extrapolation beyond the tested
scope.

**Phrases to adapt:**

* "…has been incompletely characterized and is infrequently accounted for in the internal validation of…"
* "The effect size is small for the majority of features and is absent for most features with a clear [X] basis…"
* "…external validation of every finding can be impractical/infeasible."

---

### 3.2 Venet, Dumont & Detours 2011 — the random-control battery

**Citation.** Didier Venet, Jacques E. Dumont, Vincent Detours. "Most random gene expression
signatures are significantly associated with breast cancer outcome." *PLoS Computational Biology*
**7**(10): e1002240 (2011).
DOI: <https://doi.org/10.1371/journal.pcbi.1002240>
Full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3197658/>

> **VERIFIED** via Crossref (`/works/10.1371/journal.pcbi.1002240`) and Europe PMC (PMC3197658,
> `isOpenAccess: Y`, licence CC BY). Confirmed from source: title, authors **Venet D, Dumont JE,
> Detours V**, year, journal, volume 7, article e1002240, DOI, PMID 22028643.
> Received 2011-04-27; Accepted 2011-09-07 (≈4.4 months to acceptance) — from the retrieved record.
> Section order, figure/table/reference counts and word count computed from the retrieved JATS XML.
>
> **Correction logged:** the author list was initially mis-recalled during compilation as "Venet,
> Dhahbi, Delorenzi". The retrieved record shows **Dumont** and **Detours**. The draft's
> §2.7 lists "Venet et al. 2011" as `[UNVERIFIED]` — it is now verified, with the corrected authors.

**Why it is structurally analogous.** P1's finding that ~76–82% of per-target "pathway" signal is
reproduced by covariate-matched random gene sets is, in its logical form, this paper. Venet et al.
show that random gene signatures predict breast cancer outcome nearly as well as published ones,
identify the confounder that explains it (proliferation, via a meta-PCNA metagene), supply an
adjustment for it, and then re-score the published literature through the adjustment. That is
precisely P1's arc: random control → confound identified → adjustment shipped → prior claims
re-scored. It is also the best available model for **a paper whose contribution is a negative-control
battery**, since one of its Results subsections is literally titled "Association of negative control
signatures with overall survival."

**Actual section order** (from JATS `<title>` elements):

1. Abstract
2. Author Summary
3. Introduction
4. Results
   - Most signatures not biologically related to cancer are statistically associated with breast cancer outcome
     - *Association of negative control signatures with overall survival*
     - *Most published signatures are not significantly better outcome predictors than random signatures of identical size*
   - Meta-PCNA integrates most of the outcome-related signal contained in the breast cancer transcriptome
     - *Meta-PCNA adjustment decreases the prognostic abilities of published signatures*
     - *Most prognostic transcriptional signals are correlated with meta-PCNA*
     - *Purging cell cycle genes from a signature does not rule out proliferation signals*
   - Results are reproducible across cohorts and end-points
5. Discussion
6. Methods (Software setup; Code and data availability; Expression data; Literature signatures; Meta-PCNA index; Adjusting data for the meta-PCNA index; Association of signatures with outcome)
7. Supporting Information

**Length and figures.** 6 figures, 1 table, 68 references. Body text ≈6,800 words including Methods.

**Structural lesson.** The Results are ordered *negative control first, mechanism second, adjustment
third, robustness fourth.* The paper never opens with its instrument. It opens by showing that
signatures of postprandial laughter and of mouse social defeat predict breast cancer outcome — an
absurd-on-its-face result that makes the reader want the explanation. P1 has an equivalent hook
available and is currently not using it: **chance is 0.147, not 0**, and a site confound smeared
across 256 axes is invisible to per-axis screening. Either could open the Results the way the
laughter signature opens Venet's.

**How it handles our shared rhetorical problem.**

*How it opens the negative claim:* by reframing the statistical question rather than accusing anyone.
The load-bearing sentence, and the best single sentence found in this entire search for P1's purposes:

> "Nominal p-values do not answer the appropriate statistical question: the question is not whether a
> given set of genes is related to survival, but whether it is more related to survival than random
> sets of genes."

This *is* P1's "chance is 0.147, not 0" argument, in the register that got it published. P1 should
write its own sentence to this template.

*How it states the negative claim without overreaching:* three separate containment moves, all worth
copying.

Containment of the causal claim:

> "Yet—we cannot stress this enough—we have not shown that proliferation is a core driving force
> behind breast cancer progression."

Containment of the practical claim:

> "Our study questions the biological interpretation of the prognostic value of published breast
> cancer signatures, but has no bearing on their usefulness in the clinic: a marker may be accurate
> without yielding interesting biological insight regarding the mechanism of disease progression."

Containment of the methodological claim (they concede their own analysis choice could be wrong, and
then explain why it doesn't matter):

> "The choice of association method is of course important, as there is a possibility that it misses
> some signals captured by specific combinations of signatures and models. However, most papers use
> similarly simple methods as ours."

That second quote is the template for P1's most delicate sentence. P1 is not claiming that
morphology→molecular prediction is worthless; it is claiming that the *reported effect sizes* are
measured against the wrong baseline and that the *biological interpretation* does not follow. Venet
et al. show exactly how to separate those two claims in one sentence so that reviewers cannot collapse
them.

*How it closes:* an explicitly numbered list of four claims — "in conclusion, we have shown that 1)…
2)… 3)… 4)…". For a paper with many separately-defensible findings (P1 has at least eight), an
enumerated conclusion is the correct form. It also makes each claim individually attackable, which is
a feature: reviewers who dislike claim 3 can no longer reject claims 1, 2 and 4 with it.

**Phrases to adapt:**

* "…do not answer the appropriate statistical question: the question is not whether X, but whether X is more Y than [random / chance]."
* "Yet—we cannot stress this enough—we have not shown that…"
* "Our study questions the [interpretation] of X, but has no bearing on [its practical use]."
* "Our results also imply that such markers should be evaluated against the outcome association of comparable negative control markers."

---

### 3.3 Eklund, Nichols & Knutsson 2016 — recalibrating a baseline and deflating a literature

**Citation.** Anders Eklund, Thomas E. Nichols, Hans Knutsson. "Cluster failure: Why fMRI inferences
for spatial extent have inflated false-positive rates." *Proceedings of the National Academy of
Sciences* **113**(28): 7900–7905 (2016).
DOI: <https://doi.org/10.1073/pnas.1602413113>
Record: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4948312/>

> **VERIFIED** via Crossref (`/works/10.1073/pnas.1602413113`) and Europe PMC (PMC4948312).
> Confirmed from source: title, authors, year, journal, volume 113, pages 7900–7905, DOI,
> PMID 27357684, PMCID.
> Section order, figure/table/reference counts and the quoted passages were retrieved from the PMC
> article page. Note: the Europe PMC `fullTextXML` endpoint returns 404 for this record
> (`isOpenAccess: N`); the structure and quotes below come from the PMC article page, which does
> render the full text.
> A separate Crossref record exists for a *Correction* to this paper
> (DOI `10.1073/pnas.1612033113`) — noted for accuracy; not otherwise used here.

**Why it is structurally analogous.** This is the canonical modern paper whose entire contribution is
"the baseline everyone quotes against is wrong, therefore a large published literature is inflated."
It is P1's §"chance is 0.147, not 0" argument at full scale. It is also the exemplar for the *short*
form: if P1 goes to a length-capped venue, this shows the negative claim can be carried in ~5 pages
with 2 figures.

**Actual section order** (retrieved from the PMC article page):

1. Significance
2. Abstract
3. Introduction
4. Results
   - Comparison of Empirical and Theoretical Test Statistic Distributions
   - Spatial Autocorrelation Function of the Noise
   - Spatial Distribution of False-Positive Clusters
   - Impact on a Non-Null, Task Group Analysis
   - Permutation Test for One-Sample t Test
5. Discussion
   - Why Is Clusterwise Inference More Problematic than Voxelwise?
   - Why Does AFNI's Monte Carlo Approach, Unreliant on RFT, Not Perform Better?
   - Suitability of Resting-State fMRI as Null Data for Task fMRI
   - The Future of fMRI
6. Materials and Methods

**Length and figures.** 2 figures, 1 table, 44 references. PNAS research article, 7900–7905 = 6 pages.
Main-text word count NOT VERIFIED (the page renders full text but a reliable body-only word count was
not extracted).

**Structural lesson — the Discussion is where the work is.** Note that three of the four Discussion
subsections are *pre-emptive rebuttals framed as questions*: why is this method worse than that one,
why doesn't the alternative approach fix it, and is our null data actually null. For a negative paper,
the Discussion is not a summary; it is a defence structured as a list of the objections a hostile
reviewer will raise. **P1 should restructure its Limitations section this way** — as named questions
with answers, not as a list of caveats.

**How it handles our shared rhetorical problem.**

*How it opens:* with the scale of what is at stake, quantified, in the Significance statement:

> "Functional MRI (fMRI) is 25 years old, yet surprisingly its most common statistical methods have
> not been validated using real data. Here, we used resting-state fMRI data from 499 healthy controls
> to conduct 3 million task group analyses."

And in the Introduction:

> "Since its beginning more than 20 years ago, functional magnetic resonance imaging (fMRI) has become
> a popular tool for understanding the human brain, with some 40,000 published papers according to
> PubMed. Despite the popularity of fMRI as a tool for studying brain function, the statistical
> methods used have rarely been validated using real data."

The rhetorical machinery: (i) a hard number for the size of the affected literature, (ii) "surprisingly
… have not been validated", (iii) a hard number for the scale of their own evidence. P1 has all three
available — the size of the morphology→molecular literature, the observation that its baselines have
not been validated, and its own control battery count.

*How it presents controls:* the whole paper is a negative control. They use resting-state data as null
data and measure the false-positive rate against the nominal 5%. P1's must-FAIL battery is the same
instrument.

*How it words its limitation — and this is the most transferable sentence in the paper:*

> "One possible criticism of our work is that resting-state fMRI data do not truly comprise null data,
> as they may be affected by consistent trends or transients, for example, at the start of the
> session."

The construction "One possible criticism of our work is that X" states the objection in its strongest
form, in the authors' own voice, before a reviewer can. **P1 has one control that failed as it should
not have.** That failure should be introduced with exactly this construction. Burying it will cost far
more than surfacing it — and the draft already appears to surface it (§4.1 "The instrument's own
failure, and why it is the motivating example"), which is the right instinct.

*How it closes:* forward-looking and constructive rather than triumphant —

> "Such shared data provide enormous opportunities for methodologists, but also the ability to revisit
> results when methods improve years later."

**Phrases to adapt:**

* "…is N years old, yet surprisingly its most common [methods] have not been validated using real data."
* "One possible criticism of our work is that…"
* "…but also the ability to revisit results when methods improve."

---

### 3.4 Supporting citations for P1 (verified, but not offered as structural exemplars)

These are verified and useful to cite; their full text was not retrieved, so no claim about their
section order or internal structure is made here.

* **DeGrave, Janizek & Lee 2021** — "AI for radiographic COVID-19 detection selects shortcuts over
  signal." *Nature Machine Intelligence* **3**, 610–619.
  DOI: <https://doi.org/10.1038/s42256-021-00338-7>.
  *VERIFIED via Crossref: title, three authors, journal, volume, pages, year, DOI. Full text NOT
  RETRIEVED (not in Europe PMC).* Use as the canonical shortcut-learning precedent and as evidence
  that Nature Machine Intelligence publishes purely corrective work.

* **Roberts et al. 2021** — "Common pitfalls and recommendations for using machine learning to detect
  and prognosticate for COVID-19 using chest radiographs and CT scans." *Nature Machine Intelligence*
  **3**, 199–217.
  DOI: <https://doi.org/10.1038/s42256-021-00307-0>.
  *VERIFIED via Crossref title query and Europe PMC (author list Roberts M, Driggs D, Thorpe M,
  Gilbey J, Yeung M, Ursprung S, Aviles-Rivero AI, Etmann C, McCague C, Beer L, Weir-McCall JR,
  Teng Z, Gkrania-Klotsas E, Rudd JHF, Sala E, Schönlieb C). Full text NOT RETRIEVED.*
  **Correction logged:** a DOI ending `-00307-1` was tried first and does not resolve; the correct DOI
  is `-00307-0`. Do not cite this paper from memory.

* **Snoek, Miletić & Scholte 2019** — "How to control for confounds in decoding analyses of
  neuroimaging data." *NeuroImage* **184**, 741–760.
  DOI: <https://doi.org/10.1016/j.neuroimage.2018.09.074>. PMID 30268846.
  *VERIFIED via Crossref and Europe PMC: title, three authors, journal, volume, pages, year, DOI,
  PMID, licence recorded as CC BY. Full text NOT RETRIEVED — ScienceDirect and bioRxiv both returned
  403 to automated retrieval, so no section order is claimed.* This is the neuroimaging
  confound-regression critique the brief asked about; it is a genuine near-relative of CALIBRA (it
  evaluates existing confound-control procedures and shows where they fail). **Recommend a manual
  read before drafting** — it may well displace Eklund as exemplar 3 for the instrument-design
  sections, but that judgement cannot be made from metadata alone.

* **Jiang et al. 2011** — "Synthetic spike-in standards for RNA-seq experiments." *Genome Research*
  **21**, 1543–1551. DOI: <https://doi.org/10.1101/gr.121095.111>. PMID 21816910.
  *VERIFIED via Crossref and Europe PMC: title, eight authors, journal, volume, pages, year, DOI,
  PMID. Full text NOT RETRIEVED (Europe PMC `fullTextXML` 404).* Cite as the genomics precedent for
  injection/spike-in certification of a measurement floor.

* **Lipsitch, Tchetgen Tchetgen & Cohen 2010** — "Negative controls: a tool for detecting confounding
  and bias in observational studies." *Epidemiology* **21**, 383–388.
  DOI: <https://doi.org/10.1097/ede.0b013e3181d61eeb>. PMID 20335814.
  *VERIFIED via Crossref and Europe PMC: title, three authors, journal, volume, pages, year, DOI,
  PMID. Full text NOT RETRIEVED.* Cite as the conceptual source for a negative-control battery. Note
  it is a short commentary, so it is a citation, not a structural model.

* **Leek et al. 2010** — "Tackling the widespread and critical impact of batch effects in
  high-throughput data." *Nature Reviews Genetics* **11**, 733–739.
  DOI: <https://doi.org/10.1038/nrg2825>. *VERIFIED via Crossref: title, nine authors, journal,
  volume, pages, year, DOI.*

### 3.5 Bonus — three citations the draft flags as `[UNVERIFIED]` are now verified

From `paper/P1_CALIBRA_DRAFT.md` §2.7:

* **Venet et al. 2011** → verified, **with corrected authors** (Venet D, **Dumont JE, Detours V**).
  See §3.2.
* **Leek & Storey 2007** → "Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable
  Analysis", *PLoS Genetics* **3**, e161 (2007), DOI <https://doi.org/10.1371/journal.pgen.0030161>.
  *VERIFIED via Crossref.*
* **Johnson et al. 2007** (ComBat) → W. Evan Johnson, Cheng Li, Ariel Rabinovic, "Adjusting batch
  effects in microarray expression data using empirical Bayes methods", *Biostatistics* **8**,
  118–127, DOI <https://doi.org/10.1093/biostatistics/kxj037>. *VERIFIED via Crossref. Note: Crossref
  records the issue date as 2006 (online-first); the journal issue is 2007. Check which the target
  venue's style expects.*
* **Kather et al. 2020** → "Pan-cancer image-based detection of clinically actionable genetic
  alterations", *Nature Cancer* **1**, 789–799 (2020),
  DOI <https://doi.org/10.1038/s43018-020-0087-6>. *VERIFIED via Crossref.* If the intended reference
  was instead the MSI paper: Kather et al., "Deep learning can predict microsatellite instability
  directly from histology in gastrointestinal cancer", *Nature Medicine* **25**, 1054–1056 (2019),
  DOI <https://doi.org/10.1038/s41591-019-0462-y>. *Also VERIFIED via Crossref.*

Still `[UNVERIFIED]` and not checked in this pass: Kömen et al. 2024, Carloni et al. 2025,
Schmitt et al. 2021, Muirhead 1982 / Anderson 2003 chapter and page.

**Also resolves the `[CITATION NEEDED]` in draft §2.6** — the effective-rank-as-quality-proxy
citation. Use both:
* Olivier Roy, Martin Vetterli, "The effective rank: A measure of effective dimensionality",
  *EUSIPCO* 2007, pp. 606–610. <https://ieeexplore.ieee.org/document/7098875/>
  *VERIFIED via DBLP: title, both authors, year, venue, pages.*
* Garrido, Balestriero, Najman & LeCun, RankMe — see §5.0.

---

## 4. P2 — venues

**Venue class verdict: journal, not conference — specifically TMLR.** The reason is not prestige or
fit; it is that TMLR is the only venue surveyed whose written acceptance criteria *forbid* the exact
grounds on which a correct negative result normally dies. Conferences have improved sharply in 2026
(NeurIPS now has a formal Negative Results contribution type), but they have improved by *raising the
bar* for negative results, not by protecting them.

### Rank 1 — **TMLR** (Transactions on Machine Learning Research)

**Scope fit.** Direct. TMLR's stated scope includes "accounts of applications of existing techniques
that shed light on the strengths and weaknesses of the methods" — P2 is literally that.
<https://jmlr.org/tmlr/editorial-policies.html>

**Negative results — the decisive evidence.** TMLR has no clause naming "negative results". It has
something better: its criteria delete the vocabulary used to reject them. Verbatim from
<https://jmlr.org/tmlr/acceptance-criteria.html>:

> "Are the claims made in the submission supported by accurate, convincing and clear evidence? This
> is the most important criterion."

> "Crucially, it should **not** be used as a reason to reject work that isn't considered
> 'significant' or 'impactful' because it isn't achieving a new state-of-the-art on some benchmark.
> Nor should it form the basis for rejecting work on a method considered not 'novel enough', as
> novelty of the studied method is not a necessary criteria for acceptance. We explicitly avoid these
> terms ('significant', 'impactful', 'novel'), and focus instead on the notion of 'interest'."

And from <https://jmlr.org/tmlr/editorial-policies.html>:

> "Papers should be accepted if they meet the criteria, even if the contribution or significance of
> the work is modest."

**Published precedent at this venue.** Gurukar et al., "Benchmarking and Analyzing Unsupervised
Network Representation Learning and the Illusion of Progress", TMLR 2022,
<https://openreview.net/forum?id=GvF9ktXI1V>. Also Tönshoff, Ritzert, Rosenbluth & Grohe, "Where Did
the Gap Go? Reassessing the Long-Range Graph Benchmark", TMLR 2024,
<https://openreview.net/forum?id=Nm0WX86sKv>.

**Format.** "Submissions may be any length, but a paper's length should be justified by its content."
TMLR style file mandatory; appendix after references; double-blind; arXiv preprints explicitly
permitted provided the submission does not link to the named version.

**OA / APC.** Free both ways: "TMLR imposes no fees or payments to authors, reviewers, action editors,
or editors-in-chief." CC BY 4.0, authors retain copyright.

**Time to first decision.** Target ≈**9 weeks** end-to-end (reviews due in 2 weeks for ≤12-page main
bodies, 4 weeks above; discussion 2–4 weeks). <https://www.jmlr.org/tmlr/faq.html>

**Bonus route to conference visibility.** TMLR papers awarded a J2C / Featured / Outstanding
certification are eligible for the joint NeurIPS/ICLR/ICML **Journal-to-Conference Track** (150-paper
capacity per conference). <https://neurips.cc/public/JournalToConference> This gets P2 conference
exposure without conference reviewing.

**Flag.** TMLR does not accept expanded versions of papers already published with proceedings. If P2
goes to a workshop with archival proceedings first, TMLR is closed. Non-archival workshops and arXiv
are fine.

### Rank 2 — **NeurIPS Evaluations & Datasets track** (formerly Datasets & Benchmarks)

**Scope fit.** Near-literal, and the track was renamed in 2026 precisely to make evaluation "an object
of scientific study in its own right". Effective rank is a label-free evaluation proxy; P2 is a
metric-validity paper.

**Negative results — the most explicit "yes" found anywhere.** Verbatim from
<https://neurips.cc/Conferences/2026/EvaluationsDatasetsFAQ>:

> "My paper highlights a negative result. Is ED a suitable track? Negative results, as long as they
> bring new insights and are thoroughly demonstrated via empirical evaluations, are welcome in ED
> track."

The CFP scope list includes "Present negative results, critical analyses, and use-case-inspired
evaluations" and states "Submissions need not introduce a new model or outperform prior work."
<https://neurips.cc/Conferences/2026/CallForEvaluationsDatasets> The reviewer guidelines define a
contribution type — "Evaluation Methodology and Metrics" — that is P2's exact category, with:
"Originality may be achieved by novel metrics or by showing how different evaluation assumptions lead
to different scientific conclusions. New framing is sufficient — no need to beat a baseline."
<https://neurips.cc/Conferences/2026/EvaluationsDatasetsReviewerGuidelines>

**Format.** Nine content pages; references, appendices and checklist excluded; double-blind by default
from 2026; arXiv preprints do not cause rejection; code release not mandatory for analytical papers.

**Timing — this is why it is rank 2 and not rank 1.** The NeurIPS 2026 cycle closed in May 2026
(<https://neurips.cc/Conferences/2026/Dates>). NeurIPS 2027 dates are **NOT VERIFIED** (not yet
published). Choosing this means waiting roughly nine months to submit.

### Rank 3 — **ICLR 2027** — the only conference deadline actually open

**Scope fit.** Strong: the CFP subject areas explicitly list self-supervised representation learning
and "visualization or interpretation of learned representations", and both RankMe's successor (LiDAR,
ICLR 2024) and the closest published precedent to P2 (Fang et al. 2024, see §5.0) are ICLR papers.

**Negative results.** No explicit policy. Reviewer guidance is sympathetic but unenforceable:
"not all work needs to compete on an established leaderboard … ask yourself whether the results are
likely to be interesting or surprising to a broad audience"; "this does not necessarily require
state-of-the-art results". <https://iclr.cc/Conferences/2027/ReviewerGuidelines>

**Format.** 9 pages at submission (10 for camera-ready), strictly enforced — over-length is a desk
reject. Unlimited references and appendices. Double-blind; arXiv permitted.

**Dates (retrieved from the official site).** Abstract **18 Sep 2026 AOE**; full paper **25 Sep 2026
AOE**; reviews 5 Nov 2026; decisions **16 Dec 2026**; conference 26–28 Apr 2027.
<https://iclr.cc/Conferences/2027/Dates>

**Flags.** (i) Reciprocal reviewing: at least one author must register to review ≥3 papers and hold a
prior publication at a listed venue, or you are capped at one submission. (ii) **Withdrawn and
rejected submissions remain permanently public and are de-anonymised on withdrawal** — a real
consideration for a paper that names a widely used metric as invalid.

### Rank 4 — **NeurIPS main track**, Contribution Type: *Negative Results*

Since 2026 NeurIPS asks authors to declare a contribution type, one of which is, verbatim:

> "**Negative Results:** The main contribution is in understanding a negative result. The bar for
> these submissions is expected to be high."
> <https://neurips.cc/Conferences/2026/MainTrackHandbook>

The reviewing rubric is genuinely well designed for P2 and worth reading before drafting — it
requires the result be "grounded in deeper analysis" rather than "simply an empirical observation that
some experiment did not turn out as expected", requires it be "surprising or unexpected … run counter
to a popularly held understanding", and explicitly states "A Negative Results paper need not identify
a path to mitigate the negative finding". <https://neurips.cc/Conferences/2026/ReviewerGuidelines>

**The catch, and it is stated twice:** "The significance bar for Negative Results papers is high" and
"The originality bar for Negative Results papers is high." NeurIPS's own FAQ implies P2 belongs in
E&D rather than main track, because its contribution is empirical evaluation rather than theory.
2027 dates NOT VERIFIED.

**One clause here is worth quoting in P2's rebuttals wherever it submits.** The Use-Inspired
guidelines instruct reviewers: "Expect 'non-standard' datasets. Submissions in this Contribution Type
will frequently be evaluated on real-world datasets that fall outside common ML benchmarks. Such
datasets should be encouraged if they are justified in relation to the use case." That is a
venue-published pre-emption of the "you didn't run this on ImageNet" objection.

### Rank 5 — **DMLR** (Journal of Data-centric Machine Learning Research) — *not the TMLR substitute it looks like*

Free, no page limit, JMLR-family process — superficially a TMLR clone. Three reasons it is worse for
P2: (i) **scope mismatch**, DMLR is about data, datasets and benchmarks, not representation metrics
(<https://data.mlr.press/submissions.html>); (ii) it **reinstates the exact word TMLR deletes** —
"the contribution must meet **a significance bar** with respect to its scientific value"
(<https://data.mlr.press/acceptance-criteria>); (iii) review is **single-blind**, so author identity
is disclosed. Reviews are requested in four weeks with no stated end-to-end target. Use only if TMLR
declines.

### Rank 6 — **ICML** — highest structural risk of the conferences

Topically fine (its own topic list includes "evaluation (methodology, meta studies, replicability and
validity…)", and RankMe was published here, which has a certain appeal). But ICML is the venue where
the incentive against P2 is *mechanised*: the reviewer form scores **Originality and Significance
separately on 1–4 scales**, and the accept/reject scale is written in terms of "high impact on at
least one sub-area of AI" (<https://icml.cc/Conferences/2026/ReviewerInstructions>). The CFP puts
novelty in the entry condition: "original and rigorous research of significant interest"
(<https://icml.cc/Conferences/2026/CallForPapers>). 8-page hard limit; over-length is automatic
rejection. 2027 dates NOT VERIFIED.

**This is the cleanest documentary contrast available:** ICML requires reviewers to put a number on
novelty; TMLR forbids reviewers from using novelty as grounds for rejection at all.

**Not the Position track.** ICML Position papers must "make an argument for a viewpoint or perspective
about what *should* be done, in contrast to main track papers, which report on advances that have
already been accomplished", and "Papers that describe technical research without advocating a real
position … are not responsive to this call".
<https://icml.cc/Conferences/2026/CallForPositionPapers> P2 reports completed results; it is a
research paper. NeurIPS's position track says the same more bluntly.

### Venues to avoid, with reasons

* **JMLR — disqualifying scope rule, do not submit.** Verbatim from
  <https://www.jmlr.org/author-info.html>: *"JMLR publishes papers on the theory and methods of
  machine learning but **does not publish applications of machine learning to other domains**. JMLR
  favors papers of interest to a broader machine learning audience and may deem a paper unsuitable if
  the editorial board finds its audience too narrow."* P2's entire evidentiary base is a
  histology/transcriptomics application. This is the one venue with a written rule that disqualifies
  P2 on precisely the ground that makes it interesting. See §1.3.

* **ReScience C — structurally ineligible.** It genuinely welcomes negative results (its `¬Re` / `¬Rp`
  title prefixes encode them). But it publishes independent *reimplementations of someone else's
  published study*, and states plainly: "Can I submit the replication of my own research? No."
  <https://rescience.github.io/faq/>

* **ICBINB ("I Can't Believe It's Not Better") — right ethos, wrong instrument.** Its mission
  statement is the most pro-negative-result prose of any venue surveyed ("re-value unexpected negative
  results", <http://icbinb.cc/>), but: no open call at time of writing, the current edition's theme is
  LLM-specific, the limit is **4 pages**, and archival status is explicitly unresolved. Four pages
  cannot carry four dissociations plus a seed-stability analysis, and archival publication there would
  **disqualify P2 from TMLR**. If used at all, use the non-archival track.

* **Nature Machine Intelligence — format and cost.** Articles are capped at **3,500 words** of main
  text with ~6 display items; the Gold OA APC is **£9,390 / $12,850 / €10,850**
  (<https://www.nature.com/natmachintell/submission-guidelines/publishing-options>). No negative-results
  policy found. Review timeline NOT VERIFIED.

* **Patterns (Cell Press) — plausible, costly, poorly documented.** APC **$4,900** per the DOAJ record
  for ISSN 2666-3899 (not verified against Cell Press's own page; cell.com blocks automated
  retrieval). Its one genuine advantage is cross-domain framing, and it has an on-genre precedent:
  Kapoor & Narayanan, "Leakage and the reproducibility crisis in machine-learning-based science",
  *Patterns* 4(9):100804 (2023), DOI <https://doi.org/10.1016/j.patter.2023.100804> (verified via
  Europe PMC: PMID 37720327, PMCID PMC10499856). Format and review timeline NOT VERIFIED.

### Recommended sequence for P2

1. **Submit to TMLR now.** Rolling, free, no length limit, double-blind, arXiv-compatible, ~9 weeks to
   decision, and the only venue that contractually forbids the "not novel enough" rejection.
2. **If certified, request the Journal-to-Conference Track** for NeurIPS/ICLR/ICML presentation.
3. **Alternative if a conference deadline is preferred this year:** ICLR 2027, abstract 18 Sep 2026.
   Note this is *either/or* with TMLR — TMLR forbids parallel archival submission.
4. **Best pure-fit conference is NeurIPS Evaluations & Datasets**, but that means waiting for the 2027
   cycle.

---

## 5. P2 — structural exemplars

### 5.0 First, the antagonist — verify this citation and get it exactly right

P2 exists in opposition to this paper, so its citation must be flawless.

**Quentin Garrido, Randall Balestriero, Laurent Najman, Yann LeCun. "RankMe: Assessing the Downstream
Performance of Pretrained Self-Supervised Representations by Their Rank." *Proceedings of the 40th
International Conference on Machine Learning* (ICML 2023), PMLR 202, pp. 10929–10974.**
<https://proceedings.mlr.press/v202/garrido23a.html> · arXiv:2210.02885

> **VERIFIED** via DBLP (`ICML` 2023, pages 10929–10974, PMLR URL) and the arXiv API
> (arXiv:2210.02885v3, `journal_ref: The Fortieth International Conference on Machine Learning, 2023,
> Honolulu, United States`). Confirmed: title, all four authors and order, year, venue, page range.

Related and also verified, because P2 must situate itself against it:

**Vimal Thilak, Chen Huang, Omid Saremi, Laurent Dinh, Hanlin Goh, Preetum Nakkiran, Joshua M.
Susskind, Etai Littwin. "LiDAR: Sensing Linear Probing Performance in Joint Embedding SSL
Architectures." ICLR 2024.** <https://openreview.net/forum?id=f3g5XpL9Kb> · arXiv:2312.04000

> **VERIFIED** via DBLP (venue `ICLR`, year 2024, OpenReview forum id `f3g5XpL9Kb`) and the arXiv API
> (arXiv:2312.04000v1, comment "Technical report"). Confirmed: title, all eight authors and order,
> year, venue. *Note:* the arXiv version is labelled a technical report and carries no journal_ref;
> the ICLR 2024 venue attribution comes from DBLP, not from arXiv.

This is a paper that proposes a *replacement* for RankMe-style rank measures on the grounds that they
are inadequate. P2 must engage with it explicitly — it is the closest thing to prior art, and a
reviewer will ask. P2's position is stronger than LiDAR's (LiDAR argues rank is a *weak* proxy; P2
argues rank and information are *dissociated*, including a case where rank is pinned at maximum while
the representation collapses to a point), but that distinction must be made in the paper, not left to
the reader.

**And the nearest neighbour of all, which P2 must not be caught unaware of:**

**Xianghong Fang, Jian Li, Qiang Sun, Benyou Wang. "Rethinking the Uniformity Metric in
Self-Supervised Learning." ICLR 2024.** <https://openreview.net/forum?id=3pf2hEdu8B> ·
arXiv:2403.00642

> **VERIFIED** via DBLP (venue `ICLR`, year 2024, OpenReview forum id `3pf2hEdu8B`, arXiv
> 2403.00642). Confirmed: title, all four authors and order, year, venue. Full text NOT RETRIEVED —
> the abstract was not fetched, so the summary below is drawn from a secondary report and should be
> **checked against the paper before citing**.

This is the same move P2 is making, one metric over: a widely used label-free representation-quality
measure (Wang & Isola's uniformity) is shown to fail properties it is assumed to have, including
mishandling dimensional collapse. Three consequences for P2:

1. **It is proof the genre places at a top conference** — useful when deciding between §4's options.
2. **It is the paper a reviewer will say P2 duplicates.** The differentiation is available and clean:
   Fang et al. critique uniformity on *axiomatic/property* grounds; P2 critiques effective rank on
   *empirical dissociation* grounds, with a measured information channel as the reference standard.
   P2's evidence is that the metric moves while the measured quantity does not — a different and
   arguably harder-to-dismiss form of argument. Say so explicitly in Related Work.
3. **Read it before drafting.** Because it was not fetched during this compilation, its actual
   argument structure is unverified here, and it is the one paper in this document whose content most
   directly bears on P2's novelty claim.

---

### 5.1 Leavitt & Morcos 2021 — "this widely used measure is neither necessary nor sufficient"

**Citation.** Matthew L. Leavitt, Ari S. Morcos. "Selectivity considered harmful: evaluating the
causal impact of class selectivity in DNNs." *International Conference on Learning Representations*
(ICLR) 2021, Poster. <https://openreview.net/forum?id=8nl0k08uMi> · arXiv:2003.01262

> **VERIFIED** via the OpenReview API (`/notes/search`), which returns
> `venue: "ICLR 2021 Poster"`, `venueid: "ICLR.cc/2021/Conference"`, `forum: 8nl0k08uMi`; and via the
> arXiv API (arXiv:2003.01262v3, two authors, Facebook AI Research).
> Section order, figure count and word count computed from the ar5iv rendering
> (<https://ar5iv.labs.arxiv.org/html/2003.01262>).
> Note: the arXiv record carries **no** journal_ref or venue comment — the ICLR 2021 attribution comes
> solely from OpenReview.

**Why it is structurally analogous.** This is the closest structural match P2 has, and it is not
close. A metric (class selectivity) is in wide use as an interpretive proxy for what a network is
doing. The authors do not argue about definitions; they build a knob that moves the metric directly,
move it in both directions, and measure whether the thing it supposedly indexes moves with it. It
doesn't. The conclusion is stated as a *logical* result — neither necessary nor sufficient — rather
than as a performance claim. P2 does the same thing with effective rank across four settings, and has
the stronger version of the "not sufficient" case (rank pinned at maximum, representation collapsed to
a point, retrieval below chance).

**Actual section order** (from the ar5iv rendering):

1. Abstract
2. 1 Introduction
3. 2 Related work — 2.1 Selectivity in deep learning · 2.2 Selectivity in neuroscience
4. 3 Approach — 3.1 Models and datasets · 3.2 Defining class selectivity · **3.3 A single knob to control class selectivity**
5. 4 Results — 4.1 Test accuracy is improved or unaffected by reducing class selectivity · 4.2 Does selectivity shift to a different basis set? · 4.3 Increased class selectivity considered harmful
6. 5 Discussion
7. Appendix A.1–A.15 (15 subsections)

**Length and figures.** Main text ≈5,800 words, **4 figures**, 0 tables — then a 15-subsection
appendix that is roughly as long again. Total document ≈13,600 words. This is the correct shape for
P2: a tight main text carrying the four dissociations, with the seed-stability sweep, per-setting
detail and robustness checks pushed to appendices.

**Structural lesson — §4.2 is the section P2 currently lacks.** "Does selectivity shift to a different
basis set?" is a *defeater check*: the most obvious way their result could be an artifact is that the
regulariser merely rotates selectivity off unit-aligned axes rather than removing it, so they test
that and rule it out. P2 has an exactly parallel obvious objection — *"your rank measurement is just
badly estimated / your effective-rank definition is the wrong one."* Given P2's own finding that
effective rank is seed-unstable from 9.1 to 34.1, that objection will be the first one raised. P2
needs its own §4.2, placed in the same position, ruling out "the metric is fine, your estimator
isn't."

**How it handles our shared rhetorical problem.**

*How it opens:* by taking the metric's appeal seriously before dismantling it. Verbatim:

> "This focus on individual neurons makes intuitive sense, as the tractable, semantic nature of
> selectivity is extremely alluring; some measure of selectivity in individual units is often provided
> as an explanation of 'what' a network is 'doing'."

and then the turn:

> "Finding intuitive ways of representing the workings of DNNs is essential for making them
> understandable and accountable, but we must ensure that our approaches are based on meaningful
> properties of the system."

That is the exact opening move P2 needs. Effective rank is *appealing* — it is label-free, cheap, has
a clean spectral definition, and gives a number when you have no labels. P2 should say so, in those
terms, before showing it doesn't work. A negative paper that opens by sneering at the metric loses the
reviewers who use it.

*How it justifies the method:* by naming the limitation of the standard approach and explaining why a
new instrument was needed —

> "But single unit ablation in trained networks has two critical limitations: it cannot address whether
> the presence of selectivity is beneficial, nor whether networks need to learn selectivity to function
> properly."

*How it states the negative claim without overreaching:* the claim is put in necessary/sufficient
terms, which are precise and hard to argue with, and the strongest form is hedged with "sometimes":

> "class selectivity in individual units is neither necessary nor sufficient for—and can sometimes even
> constrain—CNN performance."

P2's headline claim should be built the same way. Not "effective rank is wrong" but something of the
form *"effective rank is neither necessary nor sufficient for — and can move opposite to — measured
information content."*

*How it words its limitations:* one sentence, first thing in the relevant paragraph, conceding the
scope honestly and then converting it into future work:

> "One caveat to our results is that they are limited to CNNs trained to perform image classification.
> It's possible that our findings are due to idiosyncracies of benchmark datasets, and wouldn't
> generalize to more naturalistic datasets and tasks."

**This is the model for P2's biggest exposure.** P2's evidence comes from four settings in one
biomedical pipeline, not from a broad benchmark sweep — a reviewer will say "you showed this on your
own application." Leavitt & Morcos have the mirror-image version of that exposure (only CNNs, only
image classification) and they handle it in two sentences, without defensiveness, in the Discussion.
Note also that they turn the scope limit into an argument *for* the paper: because the metric is used
everywhere, showing it fails anywhere is informative.

*How it closes:* by generalising from the specific metric to the practice —

> "Our results make a broader point about the potential pitfalls of focusing on the properties of
> single units when trying to understand DNNs… While we consider it essential to find tractable,
> intuitive approaches for understanding complex systems, it's critical to empirically verify that
> these approaches actually reflect functionally relevant properties of the system being examined."

That last clause is the thesis of P2's conclusion, with "single units" swapped for "spectral
summaries."

**Phrases to adapt:**

* "…is extremely alluring; some measure of [X] is often provided as an explanation of 'what' a network is 'doing'."
* "…but we must ensure that our approaches are based on meaningful properties of the system."
* "X is neither necessary nor sufficient for—and can sometimes even constrain—Y."
* "One caveat to our results is that they are limited to…"
* "…it's critical to empirically verify that these approaches actually reflect functionally relevant properties of the system being examined."

---

### 5.2 Kynkäänniemi, Karras, Aittala, Aila & Lehtinen 2023 — the metric is driven by something other than what you think

**Citation.** Tuomas Kynkäänniemi, Tero Karras, Miika Aittala, Timo Aila, Jaakko Lehtinen. "The Role
of ImageNet Classes in Fréchet Inception Distance." *International Conference on Learning
Representations* (ICLR) 2023. arXiv:2203.06026
Full text used: <https://ar5iv.labs.arxiv.org/html/2203.06026>

> **VERIFIED** via the arXiv API (arXiv:2203.06026v3; five authors; comment: "ICLR 2023 camera ready.
> Code: https://github.com/kynkaat/role-of-imagenet-classes-in-fid"). Confirmed: title, all five
> authors and order, venue and year (from the authors' own camera-ready comment).
> Section order and figure/table counts computed from the ar5iv rendering. OpenReview forum id NOT
> VERIFIED (not looked up).

**Why it is structurally analogous.** FID is to generative modelling what effective rank is becoming
to representation learning: the number everyone quotes, computed automatically, rarely interrogated.
This paper shows FID is disproportionately driven by a specific artifact of its ImageNet-classifier
feature space, and that one can move FID substantially *without* moving perceptual quality. That is
the same dissociation structure as P2: move the metric, hold the thing it claims to measure fixed.

**Actual section order** (from the ar5iv rendering):

1. Abstract
2. 1 Introduction
3. 2 What does FID look at in an image? (Our visualization technique · Observations from individual images · Observations from aggregates of images)
4. 3 Probing the perceptual null space in FID (Top-1 histogram matching · Matching all fringe features · Top-N histogram matching)
5. 4 Practical example: ImageNet pre-trained GANs
6. 5 Conclusions
7. Appendices A–E, mirroring the main section titles

**Length and figures.** Main text ≈5,500 words, **7 figures and 2 tables in the main text**, plus 12
further figures and 1 table in appendices A–E. Note the appendix organisation: each appendix is titled
after the main-text section it expands, one-to-one. For P2, with four dissociation settings, that is
a clean pattern — four main-text subsections, four correspondingly-titled appendices.

**Structural lesson — "probing the perceptual null space."** Section 3's framing is worth stealing
outright. They construct the set of changes that move FID but are perceptually null. P2's third
result — *rank pinned at maximum while the representation collapses to a point* — is the same object:
the null space of the metric, exhibited constructively. P2 should consider naming a section
something like "Probing the information null space of effective rank."

**How it handles our shared rhetorical problem.**

*How it opens:* it spends the first paragraph establishing that the field is healthy and productive
and that evaluation therefore matters — "Given the large number of applications and rapid development
of the models, designing evaluation metrics for benchmarking their performance is an increasingly
important topic. It is crucial to reliably rank models and pinpoint improvements caused by specific
changes in the models or training setups." — and only then narrows to FID's dominance ("FID continues
to be the primary tool for quantifying progress"). The critique is framed as *maintenance of a shared
tool*, not as an attack.

*How it states the negative claim without overreaching — and this is the best example of the move in
either paper's literature.* The Conclusions section **opens by listing what the metric is good for**,
and only then names the failure regime:

> "The numerical values of FID have a number of important uses. Large values indicate training
> failures quite reliably, and FID appears highly dependable when monitoring the convergence of a
> training run. FID improvements obtained through hyperparameter sweeps or other trivial changes
> generally seem to translate to better (subjective) results, even when the distributions are well
> aligned.
>
> The caveats arise when two sufficiently different architectures and/or training setups are compared."

**P2's conclusion should be built on this skeleton.** The draft title already gestures at it — "the
one-sided use that survives" — which is exactly right. Effective rank is not useless: a collapsed
spectrum really does indicate a problem. What it does not do is track information content, and it
particularly does not do so when comparing across configurations. Kynkäänniemi et al. show how to say
that in a paragraph without either overclaiming or hedging the paper into meaninglessness.

*How it words its recommendation:* concrete, cheap, and adoptable —

> "As a partial solution, the FID improvements should at least be verified using a non-ImageNet trained
> Fréchet distance."

Note "As a partial solution … should at least". A negative paper is much easier to accept if it ends
with a minimal, specific, low-cost thing the reader can do on Monday. P2 needs its equivalent sentence:
if you report effective rank, at least also report [X].

*How it handles its own uncertainty:* it flags where its evidence is anecdotal, in a footnote, rather
than dressing it up —

> "Based on personal communication with individuals who have trained over 10,000 generative models."

*And it concedes the limits of its own remedy:*

> "This effect is difficult to quantify because the current widespread metrics (KID and
> Precision/Recall) also rely on the feature spaces of ImageNet classifiers."

**Phrases to adapt:**

* "The numerical values of X have a number of important uses. […] The caveats arise when…"
* "Particular care should be exercised when…, as it may compromise the validity of X as a quality metric."
* "As a partial solution, X improvements should at least be verified using…"

---

### 5.3 Musgrave, Belongie & Lim 2020 — the "reality check" template

**Citation.** Kevin Musgrave, Serge Belongie, Ser-Nam Lim. "A Metric Learning Reality Check."
*European Conference on Computer Vision* (ECCV) 2020, pp. 681–699.
DOI: <https://doi.org/10.1007/978-3-030-58595-2_41> · arXiv:2003.08505
Full text used: <https://ar5iv.labs.arxiv.org/html/2003.08505>

> **VERIFIED** via DBLP (venue `ECCV`, year 2020, pages 681–699, DOI 10.1007/978-3-030-58595-2_41)
> and the arXiv API (arXiv:2003.08505v3, three authors).
> Section order and figure/table counts computed from the ar5iv rendering.

**Why it is structurally analogous.** Included as the third exemplar because it solves a problem the
other two don't: **how to structure a paper whose Section 2 is a list of things the field does wrong.**
Its contribution is that a decade of reported improvement largely disappears under fair evaluation.
P2 is making a narrower claim about one metric, but shares the need to indict current practice without
indicting particular authors.

**Actual section order** (from the ar5iv rendering):

1. Abstract / Keywords
2. 1 Metric Learning Overview — 1.1 Why metric learning is important · 1.2 Embedding losses · 1.3 Classification losses · 1.4 Pair and triplet mining · 1.5 Advanced training methods · 1.6 Related work · **1.7 Contributions of this paper**
3. **2 Flaws in the existing literature** — 2.1 Unfair comparisons · 2.2 Weakness of commonly used accuracy metrics · 2.3 Training with test set feedback
4. 3 Proposed evaluation method — 3.1 Fair comparisons and reproducibility · 3.2 Informative accuracy metrics · 3.3 Hyperparameter search via cross validation
5. 4 Experiments — 4.1 Losses and datasets · **4.2 Papers versus reality**
6. 5 Conclusion

**Length and figures.** Main text ≈5,500 words, 3 figures, **7 tables**. Table-heavy rather than
figure-heavy — worth noting if P2's four dissociations are more legible as a table than as four plots.

**Structural lessons, three of them.**

1. **Section 1 is an overview of the field, not an introduction to the paper.** The authors establish
   shared ground and demonstrate command of the literature for six subsections before stating their
   own contribution in 1.7. For a paper telling a community it has been measuring wrong, this
   front-loaded generosity is what buys the right to the critique.
2. **"Flaws in the existing literature" is its own numbered section**, separate from both the
   introduction and the experiments, with one subsection per flaw. P2 has a version of this available:
   the ways effective rank is currently used and what each assumes.
3. **"Papers versus reality"** as a results subsection title. Blunt, memorable, and it names the
   comparison the reader actually wants.

**How it handles our shared rhetorical problem.**

*How it states the negative claim without overreaching:* the conclusion is a bulleted enumeration of
the flaws followed by a single calibrated sentence about magnitude —

> "We then ran experiments with these issues fixed, and found that state of the art loss functions
> perform marginally better than, and sometimes on par with, classic methods. This is in stark
> contrast with the claims made in papers, in which accuracy has risen dramatically over time."

"marginally better than, and sometimes on par with" is a precise, defensible quantifier. They do not
say the new methods are worthless.

*How it closes constructively:*

> "If proper machine learning practices are followed, and comparisons to prior work are done in a fair
> manner, the results of future metric learning papers will better reflect reality, and will be more
> likely to generalize to other high-impact areas like self-supervised learning."

*Limitations:* **NOT VERIFIED** — no dedicated limitations section appears in the retrieved section
list, and none was located in the conclusion. This is worth knowing as a negative datum: ECCV-format
papers of this kind often carry no explicit limitations section, which is a difference from the
journal exemplars in §3.

**Phrases to adapt:**

* "…perform marginally better than, and sometimes on par with, classic methods. This is in stark contrast with the claims made in papers…"
* "Papers versus reality" (as a section title)

---

### 5.4 Supporting citations for P2 (verified, offered as citations not exemplars)

* **Olivier Roy, Martin Vetterli. "The effective rank: A measure of effective dimensionality."**
  *EUSIPCO* 2007, pp. 606–610. <https://ieeexplore.ieee.org/document/7098875/>
  *VERIFIED via DBLP: title, both authors, year, venue, pages. Full text NOT RETRIEVED.*
  This is the definitional source. P2 must cite it, and should be careful to note that Roy & Vetterli
  proposed effective rank as a description of a matrix's spectrum — **not** as a proxy for information
  content. A useful and entirely fair framing for P2: the metric is not wrong; the *interpretation
  layered onto it downstream* is.

* **Li Jing, Pascal Vincent, Yann LeCun, Yuandong Tian. "Understanding Dimensional Collapse in
  Contrastive Self-supervised Learning."** ICLR 2022. arXiv:2110.09348
  *VERIFIED via the arXiv API: arXiv:2110.09348v3, four authors, `journal_ref: ICLR 2022`, comment
  "In Proceedings of the 10th International Conference on Learning Representations (ICLR) 2022".*
  Directly adjacent prior art: it establishes that spectral collapse is a real failure mode, which is
  *why* rank-based measures were adopted. P2's third dissociation — rank at maximum while the
  representation collapses to a point — is the counter-case to this literature and must be positioned
  against it explicitly, or a reviewer will read P2 as contradicting a well-established result when it
  is in fact identifying a case the established result does not cover.

* **Shane Barratt, Rishi Sharma. "A Note on the Inception Score."** arXiv:1801.01973.
  *VERIFIED via the arXiv API: two authors, comment "Proc. ICML 2018 Workshop on Theoretical
  Foundations and Applications of Deep Generative Models".* Note this is a **workshop** paper, not a
  main-conference paper — cite it accordingly. Useful precedent for the short-form version of the
  genre.

* **Maurizio Ferrari Dacrema, Paolo Cremonesi, Dietmar Jannach. "Are We Really Making Much Progress?
  A Worrying Analysis of Recent Neural Recommendation Approaches."** RecSys 2019. arXiv:1907.06902.
  *VERIFIED via the arXiv API: three authors, `journal_ref: Proceedings of the 13th ACM Conference on
  Recommender Systems (RecSys 2019)`.* The best-known reproducibility-critique paper in the genre;
  cite as evidence that the genre wins awards and attention rather than being buried.

---

## 6. Which single exemplar to model most closely

### 6.1 P1 → model on **Venet, Dumont & Detours 2011**

Howard et al. is the more obvious choice and it is the one to cite most often, because it is
domain-matched and because its external-validation passage licenses P1's largest declared limitation.
But it is not the one to *model*, because Howard et al. is a paper with essentially one finding
(site signatures exist and inflate accuracy) elaborated across five subsections, whereas P1 has
roughly eight separately-defensible findings and needs a structure that can carry them without the
paper reading as a list. Venet et al. is that structure. It opens on its most arresting negative
result rather than on its instrument; it orders Results as *negative control → mechanism → adjustment
→ robustness across cohorts*, which is exactly P1's arc from random gene sets through the confound to
residualisation and replication; it treats the negative-control battery as a first-class contribution
with its own named subsections rather than as an appendix of sanity checks; and it closes on an
enumerated list of four claims, which is the right form when reviewers need to be able to accept some
findings and dispute others independently. Most importantly, it has already solved P1's hardest
sentence twice — once in "the question is not whether a given set of genes is related to survival, but
whether it is more related to survival than random sets of genes", which is P1's chance-is-0.147
argument in publishable register, and once in "our study questions the biological interpretation …
but has no bearing on their usefulness in the clinic", which is how P1 separates "the reported effect
sizes are mismeasured" from "the field is worthless" in a single sentence so that a reviewer cannot
conflate the two. Take Venet's skeleton, Howard's citations and framing of the domain, and Eklund's
Discussion-as-pre-emptive-rebuttal pattern for the Limitations section.

*One caveat on this recommendation.* Hamdan et al. 2023 (§2.2) surfaced late in this compilation and
is closer to P1 in **subject** than any of the three exemplars — it is a paper about confound removal
failing in the way people do not expect. Its retrieved Results ordering (walk-through analysis →
mechanism → benchmark replication → candidate mechanisms → why it matters in practice) is a genuine
alternative skeleton, and it has the advantage of being recent and in a venue P1 is targeting. It is
not recommended as *the* model only because its contribution is single-threaded — one phenomenon,
explored in depth — whereas P1 must carry roughly eight findings, which is the specific problem
Venet's structure solves. Read Hamdan for the subject-matter framing and for how it words its central
negative claim; use Venet for the architecture.

### 6.2 P2 → model on **Leavitt & Morcos 2021**

Kynkäänniemi et al. has the better conclusion and P2 should copy that conclusion's shape almost
verbatim — open by conceding what the metric genuinely does well, then narrow precisely to the regime
where it fails, then give one cheap concrete remedy. But the *paper* to model is Leavitt & Morcos,
because the match runs the whole length of the document rather than the last page. Both papers take a
widely used, intuitively appealing, cheap-to-compute proxy; both refuse to argue about its definition
and instead build a way to move it directly and check whether the thing it supposedly indexes moves
too; both land on a necessary/sufficient formulation rather than a performance claim, which is the
form of negative claim that is hardest to rebut and easiest to state without overreaching; both need
a defeater section that rules out "your measurement is just bad" (their §4.2 on whether selectivity
merely rotates to a different basis is the direct analogue of the objection P2 will face given its own
seed-instability finding, and P2 currently has no equivalent section); and both carry the same
structural exposure — evidence drawn from one family of settings rather than a broad benchmark sweep —
which Leavitt & Morcos dispatch in two unapologetic sentences of Discussion that P2 can adapt almost
word for word. The shape is right too: roughly 5,800 words and four figures in the main text carrying
the core claim, with a fifteen-subsection appendix absorbing everything else, which is precisely how
P2's four dissociations plus seed sweep plus per-setting detail should be divided. Use Leavitt &
Morcos for the skeleton and the claim grammar, Kynkäänniemi for the conclusion and for the
"probing the null space" section framing, and Musgrave for the decision about tables versus figures.

---

## 7. Compilation corrections log

Kept because the brief asked for it and because it is the evidence that the verification protocol ran.

| # | What was initially assumed | What the source actually says | How caught |
|---|---|---|---|
| 1 | Venet et al. 2011 authors were "Venet, Dhahbi, Delorenzi" | Authors are **Venet D, Dumont JE, Detours V** | Europe PMC record for PMC3197658 |
| 2 | Roberts et al. 2021 DOI was `10.1038/s42256-021-00307-1` | That DOI does not resolve; correct DOI is **`10.1038/s42256-021-00307-0`** | Crossref `/works/{DOI}` returned no record; title query returned the correct one |

Both were caught before they reached this document. No citation in this file was written from memory.

A third, softer correction: **the widely repeated belief that "PLOS journals publish negative
results" does not transfer to PLOS Computational Biology.** The clause belongs to PLOS ONE; PLOS
Computational Biology's own criteria include "Originality" and "Innovation". Verified independently
at both URLs (§2.5). This one had not yet reached the drafts, but it is the kind of assumption that
would have.

### 7.1 What stayed unverified

Listed so that nothing here is mistaken for a confirmed fact.

**Bibliographic — full text not retrieved (citations verified, structure not claimed):** Snoek et al.
2019 (ScienceDirect and bioRxiv both 403); Jiang et al. 2011; Lipsitch et al. 2010; DeGrave et al.
2021; Roberts et al. 2021; Roy & Vetterli 2007; **Fang et al. 2024** — the last of these matters most,
because it is the paper nearest to P2's novelty claim (§5.0).

**Eklund et al. 2016:** main-text word count not extracted (section order, figure/table/reference
counts and quotes were retrieved).

**Musgrave et al. 2020:** no limitations section was located; recorded as a negative datum, not as
proof one does not exist.

**Hamdan et al. 2023:** year ambiguity between the Europe PMC `pubYear` field (2022) and the article
ID / PMID indexing (2023). Resolve before citing.

**Spisak 2022 (GigaScience):** reported from the venue search via Crossref; **not independently
re-verified** by me.

**Venue facts not verified:** conference APCs and registration requirements for NeurIPS/ICML/ICLR;
NeurIPS 2027 and ICML 2027 deadlines (not yet published); JMLR and Nature Machine Intelligence review
timelines; GigaScience APC and time to first decision (OUP renders prices via JavaScript; the ~$2,638
figure is a DOAJ record); Bioinformatics APC (same); eLife official first-decision median and APC page;
Patterns and Cell Reports Methods currency (USD inferred from an Elsevier corporate statement, not a
per-journal page); Genome Biology and npj Digital Medicine main-text word limits (not stated);
acceptance rates other than PLOS Comp Biol (33%) and Bioinformatics (~30%).

**Currency of archived text:** the *Medical Image Analysis* "significant methodological contribution"
quote comes from an Internet Archive snapshot dated 2024-02-24, and its metrics from one dated
2026-02-04, because sciencedirect.com blocks automated retrieval. Treat as probably-current, not
confirmed-current.

**Publisher self-inconsistency, unresolved:** Nature Communications states 8 days (about page) and 9
days (metrics page) to first decision; Genome Biology states 13 and 14; Bioinformatics' claimed 23-day
first decision is hard to reconcile with a ~174-day median received→accepted in its own article
records.

**Not assessed at all:** Modern Pathology and Journal of Pathology Informatics (every guideline page
returned 403/Cloudflare); npj Digital Medicine, Bioinformatics and Medical Image Analysis were
assessed only far enough to establish that they disqualify P1.

**A note on "time to first decision" generally:** most figures quoted in this document that were
derived from article records measure *submission → acceptance*, not time to first decision, and are
therefore an upper bound. Where a publisher-stated first-decision median exists it is labelled as
such, and several of those (7–10 days) are desk-screen numbers, not post-review numbers.

### 7.2 Method limitation

The WebSearch budget for this session was exhausted before research began, so this document is built
entirely on direct bibliographic-API and publisher retrieval rather than general web search. That
makes provenance stronger but coverage narrower — **this is a verified shortlist, not an exhaustive
survey.** In particular, §1.5's novelty finding for P2, and any claim of the form "no venue does X",
should be re-checked with full web search before submission.
