"""Regenerate every P2 display item, in the draft's own order.

Each figure is a standalone script; this runs them in the same interpreter so a
failure carries its own traceback, and it collects the PENDING placeholders so
that a measurement that has not landed is reported once, at the end, rather than
being noticed only by whoever squints at the artwork.

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 python make_all.py

Options
  --strict   turn every pending placeholder into a hard failure. This is what the
             test suite uses, and what a release build should use once the
             pending measurements have landed.
  --only F2  run one item by its identifier.

The figures read only from `figures/data/`, which is vendored into the
repository and never re-fetched here; `extract_from_box.py` is the only thing
that talks to the box.
"""
from __future__ import annotations

import argparse
import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import p2fig as P  # noqa: E402

#: (identifier, module, one-line description) in the order the draft presents them.
ITEMS = [
    ("F1", "fig_f1_envelope",
     "the measured retraining floor (n = 5) and the seven arm ratios inside it"),
    ("F2", "fig_f2_variance",
     "where rank's variance lives: arm against training seed"),
    ("F3", "fig_f3_floor",
     "the reproducibility floor is a property of the arm (n = 5)"),
    ("F4", "fig_f4_defeaters",
     "the four-way defeater check: the instability is in training"),
    ("F5", "fig_f5_verdict",
     "the verdict is under-determined: statistic, block, view"),
    ("F6", "fig_f6_necessity",
     "the necessity test, which went against us"),
    ("F7", "fig_f7_dilution",
     "the dilution dose-response and its miscalibration factor"),
    ("F8", "fig_f8_collapse",
     "the collapse boundary, and the 16/16 withdrawal"),
    ("T1", "tab_t1_selection_rule",
     "rank as a selection rule - and it is underpowered"),
]


def run(item: str, module_name: str, *, strict: bool) -> tuple[bool, str]:
    argv = sys.argv
    sys.argv = [module_name] + (["--strict"] if strict else [])
    try:
        module = importlib.import_module(module_name)
        importlib.reload(module)
        module.main()
        return True, ""
    except P.PendingMeasurement as exc:
        return False, f"PENDING (strict): {exc}"
    except Exception:                       # noqa: BLE001 - reported, not swallowed
        return False, traceback.format_exc()
    finally:
        sys.argv = argv


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict", action="store_true",
                    help="fail on any pending placeholder instead of drawing it")
    ap.add_argument("--only", metavar="ID", help="run one item (F1, F2, ... T1)")
    args = ap.parse_args()

    items = ITEMS if args.only is None else [i for i in ITEMS if i[0] == args.only]
    if not items:
        ap.error(f"unknown item {args.only!r}; known: {', '.join(i[0] for i in ITEMS)}")

    P.PENDING_LOG.clear()
    failures = []
    print(f"P2 figures -> {P.OUT}\n")
    for item, module_name, description in items:
        t0 = time.time()
        ok, message = run(item, module_name, strict=args.strict)
        mark = "ok " if ok else "FAIL"
        print(f"  [{mark}] {item:<3} {description}   ({time.time() - t0:.1f}s)")
        if not ok:
            failures.append((item, message))

    if P.PENDING_LOG:
        print("\nPENDING MEASUREMENTS drawn as placeholders "
              "(re-run with --strict to make these fail):")
        for line in P.PENDING_LOG:
            print(f"  - {line}")

    if failures:
        print("\nFAILURES:")
        for item, message in failures:
            print(f"\n===== {item} =====\n{message}")
        return 1

    print(f"\n{len(items)} display item(s) written as PDF, SVG and 300 dpi PNG.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
