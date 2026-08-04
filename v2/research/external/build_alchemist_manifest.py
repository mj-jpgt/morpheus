"""Build the ALCHEMIST paired-cohort manifest straight from the GDC API.

The cohort is the intersection of two GDC file sets, joined on ``case_id``:

  * ``data_type = "Slide Image"``           (program ALCHEMIST, all open, all SVS)
  * ``data_type = "Gene Expression Quantification"`` with workflow ``STAR - Counts``

GDC applies file filters *per file*, so the obvious single ``/cases`` query with
both ``data_type`` values returns 0 -- no file is both types.  The two case-ID
lists must be pulled separately and intersected, which is what this does.

Sampling rule, recorded here because it is part of the protocol:

  * one slide per patient, chosen as the lexicographically first ``file_name``.
    947 of the 1,106 paired patients have exactly one slide anyway; this rule is
    independent of file size and of tissue area, so it cannot select for
    "bigger" or "better" slides.
  * ``tissue_source_site`` is ``_missing`` for every ALCHEMIST case, so no site
    column is emitted.  ``diagnoses.primary_diagnosis`` *is* populated and is
    carried as the ``cancer`` column, which is the histology stratum the
    confound design will use.

Emits a CSV with the columns the dilution extractor already expects
(``file,size,patient,tss,cancer,slide``) plus ``case_id``, and a JSON manifest
carrying the SHA-256 of the CSV bytes and the per-endpoint counts, so the cohort
can be re-derived and checked.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

API = "https://api.gdc.cancer.gov"
PROGRAM = "ALCHEMIST"


def _post(endpoint: str, payload: dict, retries: int = 5) -> dict:
    body = json.dumps(payload).encode()
    request = urllib.request.Request(
        f"{API}/{endpoint}", data=body, headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                return json.loads(response.read().decode())
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def _files(filters: dict, fields: str) -> list[dict]:
    hits: list[dict] = []
    while True:
        page = _post("files", {"filters": filters, "size": 500, "from": len(hits),
                               "fields": fields, "sort": "file_id:asc"})
        batch = page["data"]["hits"]
        total = page["data"]["pagination"]["total"]
        if not batch:
            break
        hits.extend(batch)
        if len(hits) >= total:
            break
    return hits


def _program(*clauses: dict) -> dict:
    return {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.program.name", "value": [PROGRAM]}},
        *clauses]}


def _diagnoses() -> dict:
    """case_id -> primary_diagnosis, the only populated histology field on /cases."""
    out: dict[str, str] = {}
    while True:
        page = _post("cases", {
            "filters": {"op": "in", "content": {"field": "project.program.name", "value": [PROGRAM]}},
            "size": 500, "from": len(out), "sort": "case_id:asc",
            "fields": "case_id,submitter_id,diagnoses.primary_diagnosis"})
        batch = page["data"]["hits"]
        if not batch:
            break
        for case in batch:
            diagnoses = case.get("diagnoses") or [{}]
            out[case["case_id"]] = str(diagnoses[0].get("primary_diagnosis") or "not_reported")
        if len(out) >= page["data"]["pagination"]["total"]:
            break
    return out


def build(output_csv: Path) -> dict:
    slide_fields = "file_id,file_name,file_size,md5sum,cases.case_id,cases.submitter_id"
    expression_fields = "file_id,file_name,file_size,md5sum,cases.case_id,cases.submitter_id"
    slides = _files(_program({"op": "in", "content": {"field": "data_type", "value": ["Slide Image"]}}),
                    slide_fields)
    expression = _files(_program(
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}),
        expression_fields)

    by_case: dict[str, list[dict]] = {}
    submitter: dict[str, str] = {}
    for hit in slides:
        for case in hit["cases"]:
            by_case.setdefault(case["case_id"], []).append(hit)
            submitter[case["case_id"]] = case["submitter_id"]
    expression_cases: dict[str, list[dict]] = {}
    for hit in expression:
        for case in hit["cases"]:
            expression_cases.setdefault(case["case_id"], []).append(hit)

    paired = sorted(set(by_case) & set(expression_cases))
    diagnosis = _diagnoses()

    rows = []
    for case_id in paired:
        chosen = min(by_case[case_id], key=lambda hit: hit["file_name"])
        rows.append({
            "file": chosen["file_id"],
            "size": chosen["file_size"],
            "patient": submitter[case_id],
            # ALCHEMIST publishes no tissue_source_site; the column is kept so the
            # extractor's row contract is unchanged, and is explicitly NA rather
            # than a fabricated code.
            "tss": "NA",
            "cancer": diagnosis.get(case_id, "not_reported").replace(" ", "_").replace(",", ""),
            "slide": chosen["file_name"].split(".")[0],
            "case_id": case_id,
            "md5sum": chosen["md5sum"],
        })

    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0].keys()), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    payload = buffer.getvalue().encode()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_csv.write_bytes(payload)

    manifest = {
        "program": PROGRAM,
        "endpoint": API,
        "slide_files": len(slides),
        "slide_cases": len(by_case),
        "expression_files": len(expression),
        "expression_cases": len(expression_cases),
        "paired_cases": len(paired),
        "slides_selected": len(rows),
        "selection_rule": "one slide per paired case, lexicographically first file_name",
        "total_bytes_selected": sum(int(row["size"]) for row in rows),
        "total_bytes_all_slides": sum(hit["file_size"] for hit in slides),
        "csv_sha256": hashlib.sha256(payload).hexdigest(),
        "expression_file_ids_sha256": hashlib.sha256(
            "\n".join(sorted(hit["file_id"] for hit in expression)).encode()).hexdigest(),
        "cancer_histology_counts": {
            value: sum(1 for row in rows if row["cancer"] == value)
            for value in sorted({row["cancer"] for row in rows})},
    }
    output_csv.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="destination CSV")
    args = parser.parse_args()
    manifest = build(Path(args.output))
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
