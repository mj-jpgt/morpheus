"""F8(b): parse the vendored collapse logs into per-step tracks and endpoint pairs.

`paper/P2_FIGURES.md` F8 marks panel (b) `NEEDS EXTRACTION`. This script discharges it from the two
logs vendored under `collapse_tracks/`, and reports honestly which half of the panel has arrays and
which does not:

* ``diag_d.log`` -- the clean in-batch InfoNCE arm -- carries a **per-step** track at steps
  0/25/50/100/200/400. That is the ``12.88 -> 1.00 by step 50`` rank curve, recoverable as an array.
* ``collapse_diag.log`` -- the ``16/16`` instance -- carries **endpoint pairs only**. The script that
  wrote it probes before and after training, not on a schedule, so no array exists to extract. F8(b)'s
  own fallback applies: paired markers, labelled as endpoints, not interpolated.

**The statistic.** ``diag_d.py`` computes its ``eff-rank`` column inline as
``(sum s)^2 / sum s^2`` on the column-centred, ROW-L2-NORMALISED matrix -- that is
``RANK_VARIANTS["R3"]``, not ``spectral.CANONICAL``. Every value this script emits from that column is
therefore labelled ``R3`` in the output. No rank is computed here; the numbers are read from a log.

Run (thread caps are not optional on the shared box)::

    python3 -m morpheus.v2.research.rebase.p2.p2_f8b_tracks --output F8B_TRACKS.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

LOG_DIR = Path(__file__).resolve().parent / "collapse_tracks"

#: sha256 of each vendored log as it exists on `ubuntu@150.136.45.194`. Checked on every run: a log
#: that has drifted is not the log the published numbers came from.
SOURCES = {
    "diag_d.log": ("~/e0_run/d1_diag/diag_d.log",
                   "d2b49035e174db9b62dca95d3bed163cd6fc5e0aa8da41ec6af12231babb40f7"),
    "collapse_diag.log": ("~/e0_run/collapse_diag.log",
                          "bedadc2300d68e798d93c575b05eeb7e1a575a2362327411b104cdc4069601dc"),
}

#: The statistic of `diag_d.py:50-51`, named so it can never be quoted as the canonical one.
RANK_STATISTIC = "R3"
RANK_STATISTIC_DEFINITION = "(sum sigma)^2 / sum sigma^2 on the column-centred, row-L2-normalised matrix"
RANK_BLOCK = "raw wsi z_biology, L2-normalised rows (the same tensor the InfoNCE consumes)"

_STEP = re.compile(
    r"^\s*(?P<step>\d+)\s+loss\s+(?P<loss>[-\d.]+)\s+acc\s+(?P<acc>[-\d.]+)"
    r"\s+pos\s+(?P<pos>[-\d.]+)\s+worst-neg\s+(?P<worst_neg>[-\d.]+)"
    r"\s+min-margin\s+(?P<min_margin>[-\d.]+)\s+wsi-wsi\s+(?P<wsi_wsi>[-\d.]+)"
    r"\s+eff-rank\s+(?P<rank_r3>[-\d.]+)\s+std\s+(?P<std>[-\d.]+)\s*$")

_ARM = re.compile(r"^---\s*(?P<label>[A-Z])\.\s*(?P<description>.+?)\s*$")
#: The log aligns its columns with a variable number of spaces, and two of the six quantities --
#: including the ``16 -> 16`` rank pinning that IS the panel's subject -- use a single space. A
#: fixed-width separator silently drops them, which is how a panel ends up missing its own point.
_PAIR = re.compile(r"^\s*(?P<name>\S.*?)\s+(?P<before>[-\d.]+)\s*->\s*(?P<after>[-\d.]+)"
                   r"(?:\s+\((?P<note>[^)]*)\))?\s*$")


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def verify_sources(log_dir: Path = LOG_DIR) -> dict[str, str]:
    """Refuse to report from a log that is not byte-identical to the box original."""
    seen = {}
    for name, (origin, expected) in SOURCES.items():
        path = log_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"vendored log missing: {path} (from {origin})")
        found = _digest(path)
        if found != expected:
            raise ValueError(f"{name} has drifted from the box original {origin}: "
                             f"expected {expected}, found {found}")
        seen[name] = found
    return seen


def parse_per_step(text: str) -> list[dict[str, float]]:
    """Per-step rows of `diag_d.log`. Steps are as recorded -- nothing is interpolated."""
    rows: list[dict[str, float]] = []
    for line in text.splitlines():
        match = _STEP.match(line)
        if match is None:
            continue
        values = {key: float(value) for key, value in match.groupdict().items()}
        values["step"] = int(values["step"])
        rows.append(values)
    return rows


def parse_endpoint_pairs(text: str) -> dict[str, dict[str, object]]:
    """Arm -> {quantity: before/after}. `collapse_diag.log` holds no intermediate steps."""
    arms: dict[str, dict[str, object]] = {}
    current: dict[str, object] | None = None
    for line in text.splitlines():
        arm = _ARM.match(line)
        if arm is not None:
            current = {"description": arm.group("description"), "endpoints": {}}
            arms[arm.group("label")] = current
            continue
        if current is None or "->" not in line:
            continue
        pair = _PAIR.match(line)
        if pair is None:
            continue
        entry: dict[str, object] = {"before": float(pair.group("before")),
                                    "after": float(pair.group("after"))}
        if pair.group("note"):
            entry["note"] = pair.group("note").strip()
        current["endpoints"][pair.group("name").strip()] = entry
    return arms


def build(log_dir: Path = LOG_DIR) -> dict:
    digests = verify_sources(log_dir)
    per_step = parse_per_step((log_dir / "diag_d.log").read_text(encoding="utf-8"))
    arms = parse_endpoint_pairs((log_dir / "collapse_diag.log").read_text(encoding="utf-8"))
    return {
        "figure": "P2 F8(b)",
        "rank_statistic": RANK_STATISTIC,
        "rank_statistic_definition": RANK_STATISTIC_DEFINITION,
        "rank_statistic_source": "diag_d.py:50-51, computed INLINE; not spectral.CANONICAL",
        "rank_block": RANK_BLOCK,
        "source_digests": digests,
        "per_step_track": {
            "arm": "clean in-batch InfoNCE, lr 1e-3, all other loss terms zeroed",
            "source": "diag_d.log",
            "steps_are_as_recorded": True,
            "rows": per_step,
        },
        "endpoint_pairs_only": {
            "arms": arms,
            "source": "collapse_diag.log",
            "per_step_array_retained": False,
            "figure_instruction": ("draw as before/after paired markers labelled 'endpoint values as "
                                   "recorded; per-step array not retained'; do not interpolate"),
        },
    }


def main(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log-dir", type=Path, default=LOG_DIR)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    payload = build(args.log_dir)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    rows = payload["per_step_track"]["rows"]
    print(f"per-step track: {len(rows)} recorded steps, statistic {payload['rank_statistic']}")
    for row in rows:
        print(f"  step {row['step']:>4}  rank_{RANK_STATISTIC} {row['rank_r3']:6.2f}  "
              f"acc@1 {row['acc']:.3f}  pos {row['pos']:.4f}  worst-neg {row['worst_neg']:.4f}  "
              f"wsi-wsi {row['wsi_wsi']:.4f}")
    for label, arm in payload["endpoint_pairs_only"]["arms"].items():
        print(f"endpoint pairs, arm {label} ({arm['description']}): "
              f"{len(arm['endpoints'])} quantities, before/after only")
    return payload


if __name__ == "__main__":
    main()
