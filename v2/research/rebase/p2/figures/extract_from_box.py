"""Pull every box-resident datum the P2 figures need into `figures/data/`.

This script is the ONLY place a figure's data may enter the repository. It exists
so that "where did this number come from" has a scripted answer instead of a
remembered one: every file it writes is recorded in `data/MANIFEST.json` with its
box path, its size and its SHA-256, and every figure script reads from
`figures/data/` and nowhere else.

Two kinds of thing are pulled.

* **Vendored files** — copied byte for byte from the box (JSON readouts, run
  logs). These are listed in ``VENDORED`` and are fetched with a single
  ``tar`` stream, never with per-file ``scp``.
* **Extracted files** — small JSON summaries distilled on the box from artifacts
  that are too large to vendor (the 128 KB-per-run ``train_metrics.jsonl``
  files). The extraction is done by the inline programme in ``EXTRACTORS`` and
  the source path and SHA-256 of every input file are recorded inside the
  output, so the distillation is auditable without re-running it.

Requires ssh access to the box. It is NOT run by the test suite and NOT run by
`make_all.py`: the figures run from the vendored copies, offline.

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python extract_from_box.py
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

BOX = "ubuntu@150.136.45.194"
KEY = "~/.ssh/lambda_morpheus_nopass"

#: Files copied verbatim. Paths are relative to ``~`` on the box and are
#: reproduced under ``data/`` with the same layout, so the vendored tree reads
#: as a slice of the box tree.
VENDORED = [
    # The 2026-08-04 verified-workspace reproduction (P2_FIGURES.md path note).
    "ws_p2/out/P2_METRICS_D1.json",
    "ws_p2/out/P2_METRICS_D2.json",
    "ws_p2/out/P2_RANK_VARIANTS.json",
    "ws_p2/out/P2_ROBUSTNESS.json",
    "ws_p2/out/p2_run.log",
    # Canonicalised rank recomputation of every surviving instance.
    "ws_rank/RANK_RECOMPUTE.json",
    "ws_rank/RANK_RECOMPUTE_P1B.json",
    # D1 paired bootstraps. The UNSUFFIXED D1_PAIRED_BOOTSTRAP.json is
    # deliberately absent: it scores all 90 non-control targets, 50 of which are
    # programme_only's own supervision (P2_FIGURES.md F6 caption).
    "e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_STRATIFIED.json",
    "e0_run/d1_v2/D1_PAIRED_BOOTSTRAP_RANDOM_CONTROL.json",
    "e0_run/d1_v2/D1_PAIR_MANIFEST.json",
    # D2 readouts: retrained arms, and the re-export of the surviving seed-42
    # checkpoint that establishes readout determinism.
    "e0_run/d2_v3/D2_PER_ARTIFACT_READOUT.json",
    "e0_run/d2_v3/RECOVERED_SEED42_READOUT.json",
    # The controlled retraining envelope: five identical `programme_only` runs at
    # seed 42, read out by v2/research/rebase/d1_envelope_readout.py. This log is
    # the SOURCE OF RECORD for F1(a)/(b)/(d) -- the extractor below parses it
    # rather than recomputing anything, so the figure and the paper quote the
    # same bytes.
    "e0_run/d1_envelope_readout.log",
    # Collapse diagnostics. collapse_diag.log is the source of record for the
    # "16/16" instance and holds ENDPOINT PAIRS ONLY -- see F8(b).
    "e0_run/collapse_diag.log",
    "e0_run/d1_diag/diag_d.log",
    # Controlled short-horizon probe repeats, F4(c).
    "e0_run/d1_diag/probevar_m0.999_1.log",
    "e0_run/d1_diag/probevar_m0.999_2.log",
    "e0_run/d1_diag/probevar_m0.999_3.log",
    "e0_run/d1_diag/probevar_m0_1.log",
    "e0_run/d1_diag/probevar_m0_2.log",
    "e0_run/d1_diag/probevar_m0_3.log",
]

# --------------------------------------------------------------------------
# Extractors: run on the box, print JSON on stdout, written under data/extracted/
# --------------------------------------------------------------------------

#: F3 at n = 5. The three-seed version of this table is already in
#: ws_rank/RANK_RECOMPUTE.json under `instability_tripwire_step200_R3`; seeds 45
#: and 46 landed afterwards, in a different run directory, so the five-seed table
#: has to be built here. Key is `train_rank_tripwire_observed` WITH the `train_`
#: prefix: querying it without the prefix returns [] and [] reads as a confident
#: negative (NOTEBOOK_ENTRIES/operational_shared_box_rules_20260804T0730Z.md).
F3_EXTRACTOR = r"""
import json, hashlib, os
RUNS = {}
for s in (42, 43, 44):
    for arm in ("p", "f"):
        RUNS[f"d1_{arm}_seed{s}"] = f"/home/ubuntu/e0_run/d1_v2/d1_{arm}_seed{s}/train_metrics.jsonl"
for s in (45, 46):
    for arm in ("p", "f"):
        RUNS[f"d1_{arm}_seed{s}"] = f"/home/ubuntu/e0_run/d1_seeds4546/d1_{arm}_seed{s}/train_metrics.jsonl"
KEY = "train_rank_tripwire_observed"
EPOCH = 11
out = {
    "statistic": "R3",
    "statistic_definition": "centred + L2-row-normalised, order-2 Hill number of the singular values",
    "statistic_site": "v2/calibra/spectral.py RANK_VARIANTS['R3']; emitted by training.py as the in-run tripwire",
    "block": "training batch (in-run); states never saved, so [NOT RECOMPUTABLE] under R1",
    "key": KEY,
    "epoch": EPOCH,
    "global_step": 200,
    "runs": {},
}
for name, path in RUNS.items():
    blob = open(path, "rb").read()
    rows = [json.loads(l) for l in blob.decode().splitlines() if l.strip()]
    hits = [r for r in rows if KEY in r and r.get("epoch") == EPOCH]
    assert len(hits) == 1, (name, len(hits))
    out["runs"][name] = {
        "arm": "programme_only" if "_p_" in name else "programme_free",
        "seed": int(name[-2:]),
        "value": hits[0][KEY],
        "epoch": hits[0]["epoch"],
        "source": path,
        "source_sha256": hashlib.sha256(blob).hexdigest(),
        "source_bytes": len(blob),
    }
print(json.dumps(out, indent=1, sort_keys=True))
"""

#: F1's controlled retraining repeat. It may legitimately come back empty; the
#: figure script turns an empty result into a loud PENDING panel rather than into
#: a silent omission.
#:
#: NOTHING IS RECOMPUTED HERE. Every rank and channel value is parsed out of
#: `~/e0_run/d1_envelope_readout.log`, which was written by
#: `v2/research/rebase/d1_envelope_readout.py` importing `calibra.spectral` and
#: `calibra.residualise` -- the only place in this project a rank or channel
#: statistic may come from. The extractor records the readout module's own git
#: blob SHA-1 so the figure can assert it against the repository copy: a
#: workspace whose `spectral.py` had drifted from HEAD is exactly how the
#: variance decomposition was nearly published under a different function
#: (`NOTEBOOK_ENTRIES/WORKSPACE_DRIFT_AUDIT_ALL_20260803T2359Z.md`).
F1_REPEAT_EXTRACTOR = r"""
import json, glob, hashlib, os, re
root = "/home/ubuntu/e0_run/d1_envelope"
log_path = "/home/ubuntu/e0_run/d1_envelope_readout.log"
module = "/home/ubuntu/ws_d1/morpheus/v2/research/rebase/d1_envelope_readout.py"
out = {"root": root,
       "readout_log": log_path,
       "readout_module": "v2/research/rebase/d1_envelope_readout.py",
       "readout_module_on_box": module,
       "expected": "rep{1..5}/, rep{1..5}.npz, readout log",
       "statistic": "R1",
       "block": "rank_raw on the exported wsi_biology block; rank_residualised and "
                "channel on the cancer + pooled-TSS residualised block, top-CCA at "
                "16 components, 40 targets neither arm trained on",
       "predeclaration": "NOTEBOOK_ENTRIES/PREDECLARED_retraining_envelope_20260804T0330Z.md",
       "reps": {}, "printed_spread": {}, "complete": False}

if os.path.exists(module):
    blob = open(module, "rb").read()
    out["readout_module_blob_sha1"] = hashlib.sha1(
        b"blob " + str(len(blob)).encode() + b"\0" + blob).hexdigest()

readout = {}
if os.path.exists(log_path):
    raw = open(log_path, "rb").read()
    out["readout_log_sha256"] = hashlib.sha256(raw).hexdigest()
    text = raw.decode()
    for m in re.finditer(r"^\s*(rep\d+)\.npz\s+rank_raw=\s*([0-9.]+)\s+"
                         r"rank_resid=\s*([0-9.]+)\s+channel=([0-9.]+)\s*$", text, re.M):
        readout[m.group(1)] = {"rank_raw": float(m.group(2)),
                               "rank_residualised": float(m.group(3)),
                               "channel": float(m.group(4))}
    for label, key in (("rank (raw)", "rank_raw"),
                       ("rank (residualised)", "rank_residualised"),
                       ("channel", "channel")):
        m = re.search(r"^" + re.escape(label) + r"\s+min=([0-9.]+)\s+max=([0-9.]+)"
                      r"\s+spread=([0-9.]+)x\s*$", text, re.M)
        if m:
            out["printed_spread"][key] = {"min": float(m.group(1)),
                                          "max": float(m.group(2)),
                                          "spread": float(m.group(3))}

if os.path.isdir(root):
    for d in sorted(glob.glob(root + "/rep*")):
        if not os.path.isdir(d):
            continue
        name = os.path.basename(d)
        jl = os.path.join(d, "train_metrics.jsonl")
        rows = []
        if os.path.exists(jl):
            rows = [json.loads(l) for l in open(jl) if l.strip()]
        npz = os.path.join(root, name + ".npz")
        entry = {
            "epochs_logged": len(rows),
            "last_epoch": (rows[-1].get("epoch") if rows else None),
            "artifact_npz": npz if os.path.exists(npz) else None,
            "rank_raw": None, "rank_residualised": None, "channel": None,
        }
        if os.path.exists(npz):
            h = hashlib.sha256()
            with open(npz, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            entry["artifact_sha256"] = h.hexdigest()
        entry.update(readout.get(name, {}))
        out["reps"][name] = entry
    # Complete iff every rep has an exported artifact AND the readout reported it.
    out["complete"] = bool(out["reps"]) and all(
        r["artifact_npz"] and r["rank_residualised"] is not None and r["channel"] is not None
        for r in out["reps"].values())
print(json.dumps(out, indent=1, sort_keys=True))
"""

EXTRACTORS = {
    "extracted/F3_TRIPWIRE_STEP200_R3_n5.json": F3_EXTRACTOR,
    "extracted/F1_RETRAINING_REPEAT.json": F1_REPEAT_EXTRACTOR,
}

THREAD_CAPS = ("OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "
               "MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1")


def _ssh(command: str) -> bytes:
    return subprocess.run(
        ["ssh", "-i", KEY, BOX, command],
        check=True, stdout=subprocess.PIPE).stdout


def fetch_vendored() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    tar = _ssh("cd ~ && tar cf - " + " ".join(VENDORED))
    (DATA / "_box.tar").write_bytes(tar)
    subprocess.run(["tar", "xf", "_box.tar"], cwd=DATA, check=True)
    (DATA / "_box.tar").unlink()


def run_extractors() -> None:
    (DATA / "extracted").mkdir(parents=True, exist_ok=True)
    for rel, programme in EXTRACTORS.items():
        blob = _ssh(f"{THREAD_CAPS} python3 - <<'PYEOF'\n{programme}\nPYEOF")
        json.loads(blob)  # fail here rather than in a figure
        (DATA / rel).write_bytes(blob)


def write_manifest() -> None:
    entries = {}
    for path in sorted(DATA.rglob("*")):
        if not path.is_file() or path.name == "MANIFEST.json":
            continue
        rel = path.relative_to(DATA).as_posix()
        blob = path.read_bytes()
        entries[rel] = {
            "sha256": hashlib.sha256(blob).hexdigest(),
            "bytes": len(blob),
            "box_path": "~/" + rel if rel in VENDORED else None,
            "produced_by": None if rel in VENDORED else "extract_from_box.py",
        }
    (DATA / "MANIFEST.json").write_text(json.dumps({
        "fetched_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "box": BOX,
        "note": "Written by extract_from_box.py. Figures read from here and nowhere else.",
        "files": entries,
    }, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    fetch_vendored()
    run_extractors()
    write_manifest()
    print(f"data/ refreshed; manifest lists {len(json.loads((DATA / 'MANIFEST.json').read_text())['files'])} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
