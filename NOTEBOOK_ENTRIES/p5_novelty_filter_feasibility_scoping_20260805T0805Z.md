# P5 §4.3 — novelty-filter feasibility scoping: PubMed and Open Targets are reachable, but "no hit" is weaker evidence than it looks

**Status:** RESULT (environment scoping, not a biological measurement)
**Experiment:** P5 (planning-stage discovery engine)
**Logged:** 2026-08-05 08:05 UTC
**Predeclared:** `NOTEBOOK_ENTRIES/PREDECLARED_p5_novelty_scoping_and_pilot_funnel_20260805T0750Z.md`
**How obtained:** direct `curl` calls from this sandbox to `eutils.ncbi.nlm.nih.gov` and
`api.platform.opentargets.org` (real HTTPS round trips, no mock/replay), timed and counted with
`curl -w`. No SDK, no key.

---

## Bad news first

1. **A process deviation happened before this file existed**, and is disclosed in full in the
   predeclaration: the exact PubMed query from the task brief, one `esummary` call, a 10-request
   rate-limit burst, and one Open Targets query were run before any predeclaration was written.
   Recorded as a deviation, not hidden — see the predeclaration's §0.
2. **A "no hit" from either API is materially weaker evidence than a first read suggests**,
   especially for PubMed, and especially for the specific kind of claim this pipeline will
   surface. This is the substantive finding of this pass, not a footnote: see §3.
3. **Neither API's coverage extends to the actual claim type.** Open Targets is gene-disease
   evidence; PubMed tier-3 as specified is a generic gene-AND-cancer-AND-histology co-occurrence
   search. Neither directly indexes "does axis carry a morphology-molecular correlation for this
   gene" — the plan already says tier 3 "cannot be made airtight" (§1), and this pass sharpens
   *why*, not just restates it.

---

## What was run

Four categories of request, all from this sandbox, all real round trips:

1. The exact brief query against PubMed `esearch`:
   `("TP53"[Title/Abstract] AND "lung adenocarcinoma"[Title/Abstract] AND (histology OR morphology))`.
2. One `esummary` call for a returned PMID.
3. Two sustained-rate PubMed tests: 12 requests at ~0.34 s spacing (~2.9 req/s) and 10 requests at
   0.5 s spacing (2 req/s), plus one uncontrolled 10-request burst with no delay to find where a
   `429` appears.
4. One deliberately obscure PubMed query (`ZBTB7B AND cholangiocarcinoma AND (histology OR
   morphology)`) to see the *shape* of a genuine zero-hit response.
5. Open Targets GraphQL: one `associatedDiseases` query for TP53, then a 12-request burst with no
   delay.

## What happened

**PubMed E-utilities is reachable and usable, with a real, cleanly-signalled rate limit.**
- The brief's exact query returned HTTP 200 in 0.31 s with `count=413` and a real,
  MeSH-translated query (`histology`/`morphology` both expanded through MeSH subheadings). Not a
  mock, not an error page with a 200 status.
- `esummary` for one of the returned PMIDs returned full metadata (title, authors, journal,
  DOI, PMC ID) in 0.19 s.
- **Uncontrolled burst (0 s spacing, no key): 2 of 10 requests came back `429`** (positions 4 and
  8 of 10, in 3.05 s wall time for the batch) — this is NCBI's documented ~3 req/s anonymous
  ceiling being enforced in practice, not just in the docs.
- **At ~2.9 req/s (0.34 s spacing) and at 2 req/s (0.5 s spacing): 0 of 22 requests were
  throttled** across two separate runs (12 and 10 requests). A scriptable checker that paces
  itself at or under 3 req/s will not need a key for low-volume, first-pass use.
- **The rate limit is signalled explicitly (`429`), not silently** — no run returned a 200 with an
  empty or truncated body under load. This matters specifically because a silent throttle would be
  indistinguishable from a genuine "no hit" for an automated tier-3 checker, which is the single
  most dangerous failure mode a novelty filter could have. It was checked for and not found.
- **A genuine zero-hit query returns a clean, checkable negative**: `count="0"`, `idlist: []`, and
  an explicit `warninglist.outputmessages: ["No items found."]` — distinguishable in code from a
  malformed-query false negative (which would show a `querytranslation` that dropped a clause
  silently, not an explicit "no items" message).

**Open Targets GraphQL is reachable, fast, and returned real curated data with no observed rate
limiting at pilot volume.**
- One `associatedDiseases` query for TP53 returned in 0.94 s with five real, score-ranked
  associations (Li-Fraumeni syndrome 0.876, hepatocellular carcinoma 0.797, head and neck
  squamous cell carcinoma 0.777, choroid plexus papilloma 0.766, colorectal cancer 0.750) out of a
  reported `count: 5638` total associations for this one gene.
- A 12-request burst with **zero delay** returned 12/12 `200`s in 3.17 s (~3.8 req/s sustained,
  no throttling observed at this volume). No anonymous-key requirement is documented for this
  endpoint, consistent with what was observed.

## Technical

| check | result |
|---|---|
| PubMed `esearch`, brief's exact query | HTTP 200, 0.31 s, `count=413` |
| PubMed `esummary` | HTTP 200, 0.19 s, full metadata |
| PubMed burst, 0 s spacing, n=10 | 2/10 `429` (positions 4, 8), 3.05 s total |
| PubMed sustained, 0.34 s spacing, n=12 | 0/12 `429` |
| PubMed sustained, 0.5 s spacing, n=10 | 0/10 `429` |
| PubMed deliberate zero-hit query | `count=0`, `idlist=[]`, explicit "No items found." |
| Open Targets GraphQL, TP53 associations | HTTP 200, 0.94 s, 5 real scored rows of 5,638 |
| Open Targets burst, 0 s spacing, n=12 | 12/12 `200`, 3.17 s total, no throttling observed |

## In plain terms

Both free, keyless APIs the plan names work from this environment, respond fast, and PubMed's rate
limit shows up as an explicit, checkable error rather than a silent gap that could be mistaken for
"nothing published". A scriptable tier-2/tier-3 checker is mechanically buildable today without a
key, provided it paces PubMed calls at no more than ~3/s.

## Meaning for the claim

This licenses a narrow, specific claim: **the query mechanics for tier-2 (Open Targets) and tier-3
(PubMed) checks work from this sandbox, at low volume, without a key.** It does **not** license
"this pipeline can honestly call a finding novel." Two reasons, both already flagged in the plan
and sharpened here rather than newly discovered:

1. **PubMed tier-3 false-negative rate is a judgment call, not computed, and the honest judgment is
   that it is not small for this pipeline's actual claim type.** The brief's query form (`gene AND
   cancer-type AND (histology OR morphology OR imaging)`) is a *generic* co-occurrence search. It
   will reliably surface "is TP53 studied in lung adenocarcinoma at all" (yes, 413 hits) but is not
   shaped to find "does axis 47's morphology correlate with TP53 expression in LUAD" even when such
   a paper exists, because that is a narrow methodological claim unlikely to be captured by a
   simple title/abstract Boolean over MeSH-expanded generic terms. Add to that gene-alias mismatches
   (HGNC symbol vs. historical names), non-English literature, and unindexed preprints. A double-digit
   false-negative percentage for the pipeline's actual claim type is the right order of magnitude
   to assume, not a number this pass computed precisely — stated as the judgment call the task asked
   for, not a false precision.
2. **Open Targets tier-2 covers gene-disease associations, not gene-morphology associations.** A
   "no hit" there says "this gene is not curated as linked to this disease" — informative for tier-1/
   tier-2 target plausibility, silent on whether the specific morphology-molecular correlation this
   pipeline would surface has ever been published.

So: **usable for "is this an obvious, already-curated gene-disease link" (tier 2) and "is there any
indexed literature co-occurrence at all" (tier 3), not usable as proof of absence for the pipeline's
actual claim type.** Every hypothesis card this pipeline eventually produces needs this caveat
attached explicitly, every time, per the plan's own §1 and §6 — not summarised once here and
dropped.

## Files / commits
- This entry; predeclaration `NOTEBOOK_ENTRIES/PREDECLARED_p5_novelty_scoping_and_pilot_funnel_20260805T0750Z.md`.
- No code changed for this part; all evidence is the raw HTTP responses quoted above (not saved to
  disk as artifacts — they are timestamped, reproducible API calls, not something with its own
  provenance hash to verify against).
