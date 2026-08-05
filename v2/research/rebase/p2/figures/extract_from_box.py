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
    # The retraining floor measured for EVERY statistic, view and block the five
    # exported same-seed repeats allow -- draft 4.1's floor is canonical R1 on the
    # exported `wsi_biology` block and nothing else, and 4.1a was applying it to
    # eight other statistics and four other blocks. Produced by
    # `v2/research/rebase/p2/p2_envelope_floors.py` on a workspace verified
    # file-by-file against `git ls-tree`; the JSON carries the sha256 of each
    # `rep{n}.npz` it read and its own `absent_blocks` list, which is what makes
    # "no floor exists for this block" a recorded result rather than a silence.
    "ws_floor/out/P2_ENVELOPE_FLOORS.json",
    "ws_floor/out/floors_run.log",
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
    # The decorrelation ablation at m = 0.999, one seed per level, 400 steps.
    # F9's whole source: it prints BOTH rank statistics and the RNA-view mutual
    # cosine on the same rows, which is what makes the dissociation co-measured
    # rather than compared across tables.
    "e0_run/d1_diag/ablate_decorr0.0.log",
    "e0_run/d1_diag/ablate_decorr0.01.log",
    "e0_run/d1_diag/ablate_decorr0.04.log",
    # The momentum seed replication behind draft 5.2/5.4. Vendored so the floor
    # audit can read 10.45 and 3.18 from the bytes rather than from a notebook
    # entry. mseed_m0.999_s42 is ALSO an independent same-seed repeat of the
    # decorr = 0.04 arm of the ablation above -- same seed, same momentum, same
    # capacity, same lr, and `d1_momentum_probe.py` has no schedule that depends
    # on the step budget -- which is the only like-for-like pair this project has
    # in the momentum/probe regime. See draft 4.1a.
    "e0_run/d1_diag/mseed_m0.999_s42.log",
    "e0_run/d1_diag/mseed_m0.999_s43.log",
    "e0_run/d1_diag/mseed_m0.999_s44.log",
    "e0_run/d1_diag/mseed_m0_s42.log",
    "e0_run/d1_diag/mseed_m0_s43.log",
    "e0_run/d1_diag/mseed_m0_s44.log",
    # The six-arm learning-rate test behind draft 5.2a -- the paper's one
    # ESTABLISHED mechanism result, and until now the only §5 block whose audit
    # rows resolved against the draft's own table rather than against bytes. The
    # box filenames carry the arm's lr and momentum, which is why they are not
    # `lr_L{1..6}.log`. Each log's own first line echoes the resolved argv, so
    # vendoring them also recovers the two parameters the notebook entries never
    # recorded (decorrelation 0.04, seed 42) AND the step budget, which is the
    # column the audit reads.
    "e0_run/d1_diag/lr_L1_hi_m0.9.log",
    "e0_run/d1_diag/lr_L2_lo_m0.99.log",
    "e0_run/d1_diag/lr_L3_hi_m0.log",
    "e0_run/d1_diag/lr_L4_lo_m0.999.log",
    "e0_run/d1_diag/lr_L5_hi_m0.999.log",
    "e0_run/d1_diag/lr_L6_lo_m0.log",
    # §5.2's HEADLINE MOMENTUM SWEEP -- the step-0-to-600 table, and the source of
    # every value §5.4 row 1 and §5.2's per-step folds are read from. §6.4 flags
    # only the `lr_L*` gap, so this table's provenance was exactly as weak as
    # §5.2a's without being labelled so: audit rows 43, 46 and 47 resolved against
    # the draft's own table. The runs were launched with a 1,500-step budget and
    # the table reads them to 600.
    "e0_run/d1_diag/long_m0.log",
    "e0_run/d1_diag/long_m0.9.log",
    "e0_run/d1_diag/long_m0.99.log",
    "e0_run/d1_diag/long_m0.999.log",
    # §5.2 MEASUREMENT 2, the staleness falsification: rank AND the key-to-encoder
    # cosine on the same rows at step 100. These are a DIFFERENT run family from
    # `long_*` above -- 300-step runs, and at step 100 they read 2.58 / 6.65 / 6.89
    # where `long_*` reads 1.62 / 6.49 / 7.03 for the same nominal configuration.
    # §5.2 quotes both families and does not say so; see `known_source_disagreements`.
    "e0_run/d1_diag/mom_0_d0.04.log",
    "e0_run/d1_diag/mom_0.99_d0.04.log",
    "e0_run/d1_diag/mom_0.999_d0.04.log",
    # §5.2 MEASUREMENT 3, the capacity sweep -- the two logs draft §6.2 says are
    # "not vendored" and whose reading step it says "was never recorded". They
    # are, and it is: both arms print a row at STEP 150, where capacity 64 reads
    # 6.17 and capacity 4,096 reads 2.16, which are the two numbers §5.2 quotes.
    # A different harness (`decorr_causal.py`) with a five-column header carrying
    # a trailing decorrelation-LOSS column; `PROBE_HEADERS` names it rather than
    # guessing at it. The capacity-4,096 arm's log is the `decorr_causal_*` one
    # because that sweep held capacity at the value D1 runs and varied only
    # decorrelation.
    "e0_run/d1_diag/qsweep_d0.04_cap64.log",
    "e0_run/d1_diag/qsweep_d0.0_cap64.log",
    "e0_run/d1_diag/qsweep_d0.04_cap512.log",
    "e0_run/d1_diag/qsweep_d0.0_cap512.log",
    "e0_run/d1_diag/decorr_causal_0.04.log",
    "e0_run/d1_diag/decorr_causal_0.0.log",
    # §5.2's TURNOVER FALSIFICATION, read at step 250. Audit rows 51, 52 and 53.
    "e0_run/d1_diag/turn_cap2048_m0.9.log",
    "e0_run/d1_diag/turn_cap2048_m0.95.log",
    "e0_run/d1_diag/turn_cap4096_m0.95.log",
    "e0_run/d1_diag/turn_cap8192_m0.95.log",
    "e0_run/d1_diag/turn_cap8192_m0.99.log",
    # THE RETRAINING FLOOR ON THE FIXED HELD-OUT PROBE. Draft §6.2 predeclared
    # this run in these words -- "five same-seed repeats of the `programme_free` /
    # 500-step configuration with `d1_momentum_probe.py` attached, read at a fixed
    # step" -- because every rank number in §5 is measured on that block and no
    # floor had ever been measured for it, which made all of §5 `unjudgeable`
    # rather than passing or failing. Ten runs: five identical same-seed repeats
    # of EACH of the two arms §5 compares. Produced by `p2_probe_floors.py` on a
    # workspace verified file-by-file against `git ls-tree`; the JSON carries the
    # sha256 of every probe state it read, every per-repeat value, and its own
    # `absent` list of what the floor still does not cover.
    "e0_run/d1_probefloor/out/P2_PROBE_FLOORS.json",
    "e0_run/d1_probefloor/out/probe_floors_run.log",
]

#: The ten probe-floor logs themselves, vendored beside the JSON so the two rank
#: columns can be read from the run's own bytes and not only from the summary.
VENDORED += [f"e0_run/d1_probefloor/probefloor_m{m}_rep{r}.log"
             for m in ("0.999", "0.0") for r in range(1, 6)]

#: THE CENTRING MEASUREMENT. §6.2's last open row: the RNA-view mutual cosine
#: recomputed on the CENTRED representation for §5.2a's three `lr = 1e-3` arms,
#: which needed those runs' activations and so needed the runs re-done with
#: `export_dir` attached. Three same-seed repeats of each arm rather than the one
#: seed §5.2a has, so that "moves" and "flat" are verdicts against a measured
#: spread. The reading rule was predeclared at
#: `NOTEBOOK_ENTRIES/PREDECLARED_centred_cosine_20260804T1700Z.md` before any arm
#: ran. The JSON carries the sha256 of every state it read and the state/log
#: guard's per-row deltas.
VENDORED += ["e0_run/d1_lrcentre/out/P2_CENTRED_COSINE.json",
             "e0_run/d1_lrcentre/out/centred_cosine_run.log"]
VENDORED += [f"e0_run/d1_lrcentre/lrc_m{m}_rep{r}.log"
             for m in ("0", "0.9", "0.999") for r in range(1, 4)]

#: THE THREE FLOORS THAT CLOSE §4.1a'S THREE UNJUDGEABLE ROWS. Each was named in
#: §6.2 with what it would cost, and each is the same measurement at a setting no
#: repeat had been run at: a **step-600** budget for §5.4 row 1, a third arm at
#: **m = 0.99** for §5.4 limit 2, and **capacity 64** against capacity 4,096 for
#: §5.2 measurement 3, read at the step its own logs turn out to record. Scored
#: once per PAIR, because the two step-600 rows do not share a second arm.
VENDORED += ["e0_run/d1_probefloor600/out/P2_PROBE_FLOORS_S600_m0999_m0.json",
             "e0_run/d1_probefloor600/out/P2_PROBE_FLOORS_S600_m0999_m099.json",
             "e0_run/d1_probefloor600/out/probe_floors600_run.log",
             "e0_run/d1_capfloor/out/P2_PROBE_FLOORS_CAP.json",
             "e0_run/d1_capfloor/out/cap_floors_run.log"]
VENDORED += [f"e0_run/d1_probefloor600/pf600_m{m}_rep{r}.log"
             for m in ("0.999", "0", "0.99") for r in range(1, 6)]
VENDORED += [f"e0_run/d1_capfloor/cap_cap{c}_rep{r}.log"
             for c in ("64", "4096") for r in range(1, 6)]

#: §5.4 LIMIT 2 PUSHED UNTIL IT BROKE OR HELD. The row above closed the last
#: unjudgeable selection and closed it fragile: it clears its own floor by 5.6%
#: under R3 while the same ten runs separate the arms by only 1.138x worst case
#: under that statistic, and by 1.453x under canonical R1. It is the one selection
#: this paper makes at the hyperparameter the project actually ships, so it was
#: pushed on three axes rather than restated -- predeclared in
#: `NOTEBOOK_ENTRIES/PREDECLARED_probe_floor_n10_and_momentum_grid_20260805T0200Z.md`
#: before anything ran.
#:
#: * **repeats**, 5 -> 10 per arm at the shipped setting (`_N10`), with repeats
#:   6-10 also scored ALONE (`_LATE5`) as an independent five, because a floor
#:   that only widens because the batches differ is a batch effect wearing a
#:   floor's clothes;
#: * **statistic**, every key of `RANK_VARIANTS` beside the published
#:   alternatives, with duplicates detected numerically rather than assumed;
#: * **momentum**, a five-point grid {0, 0.98, 0.99, 0.995, 0.999} at the same
#:   step and seed (`_GRID`), so the shipped comparison is placed against its
#:   neighbours instead of being read as a two-point contrast.
#:
#: Each file also carries the same pair scored at every step the probe states
#: were exported at, which cost no GPU and is where the sharpest result is.
VENDORED += ["e0_run/d1_probefloor600/out/P2_LIMIT2_STRESS_N5.json",
             "e0_run/d1_probefloor600/out/P2_LIMIT2_STRESS_N10.json",
             "e0_run/d1_probefloor600/out/P2_LIMIT2_STRESS_LATE5.json",
             "e0_run/d1_probefloor600/out/P2_LIMIT2_STRESS_N10_RNA.json",
             "e0_run/d1_probefloor600/out/P2_MOMENTUM_GRID.json",
             "e0_run/d1_probefloor600/out/limit2_stress_run.log"]
VENDORED += [f"e0_run/d1_probefloor600/pf600_m{m}_rep{r}.log"
             for m in ("0.999", "0.99") for r in range(6, 11)]
VENDORED += [f"e0_run/d1_probefloor600/pf600_m{m}_rep{r}.log"
             for m in ("0.995", "0.98") for r in range(1, 6)]

#: THE EXPORTED-BLOCK FLOOR ON THE UNSTABLE ARM. Every exported-block floor this
#: paper quotes -- 3.295x included -- is five same-seed `programme_only` retrains,
#: the STABLE arm, which is the configuration the paper's own block-matching
#: argument says flatters a floor. `_PF` is the identical measurement on
#: `programme_free`: the same chain script with `--objective-profile` as the only
#: changed flag, same seed 42, same 40 epochs, same export. `_PO_RECHECK` is the
#: stable arm re-scored through the same invocation in the same session, so the
#: two arms are comparable by construction rather than by trust -- it reproduces
#: 3.2947x / 3.1110x / 1.0193x / 1.0200x and every per-repeat value in draft 4.1.
#: Predeclared in
#: `NOTEBOOK_ENTRIES/PREDECLARED_unstable_arm_exported_floor_20260805T0045Z.md`.
#:
#: SEEDAXIS_{f,p} are the same module on the ten exported 40-epoch D1 artifacts at
#: seeds 42-46. That is 4.2's axis, NOT a retraining floor, and it is vendored so
#: the two quantities can be told apart in one place rather than confused.
VENDORED += ["e0_run/d1_envelope_pf/out/P2_ENVELOPE_FLOORS_PF.json",
             "e0_run/d1_envelope_pf/out/P2_ENVELOPE_FLOORS_PO_RECHECK.json",
             "e0_run/d1_envelope_pf/out/d1_envelope_pf_readout.log",
             "e0_run/d1_envelope_pf/out/floors_pf_run.log",
             "e0_run/d1_envelope_pf/out/floors_po_recheck.log",
             "e0_run/pf_seedaxis/out/SEEDAXIS_f.json",
             "e0_run/pf_seedaxis/out/SEEDAXIS_p.json"]

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
        # `.gitattributes` is a repository file, not vendored evidence: this
        # directory is marked `* -text` so the digests below stay meaningful, and
        # the marker itself must not be listed as something the box produced.
        if not path.is_file() or path.name in ("MANIFEST.json", ".gitattributes"):
            continue
        rel = path.relative_to(DATA).as_posix()
        blob = path.read_bytes()
        entries[rel] = {
            "sha256": hashlib.sha256(blob).hexdigest(),
            "bytes": len(blob),
            "box_path": "~/" + rel if rel in VENDORED else None,
            "produced_by": None if rel in VENDORED else "extract_from_box.py",
        }
    # Bytes, with an explicit LF, not `write_text`: this directory is marked
    # `* -text`, so a run on Windows would otherwise rewrite every line of the
    # manifest as CRLF and present it as a change to every vendored digest.
    (DATA / "MANIFEST.json").write_bytes((json.dumps({
        "fetched_utc": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "box": BOX,
        "note": "Written by extract_from_box.py. Figures read from here and nowhere else.",
        "files": entries,
    }, indent=1, sort_keys=True) + "\n").encode("utf-8"))


def main() -> int:
    fetch_vendored()
    run_extractors()
    write_manifest()
    print(f"data/ refreshed; manifest lists {len(json.loads((DATA / 'MANIFEST.json').read_text())['files'])} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
