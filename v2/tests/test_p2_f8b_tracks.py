"""F8(b)'s extraction, pinned: the arrays that exist, the endpoints that are all that exists.

`paper/P2_FIGURES.md` F8 marks panel (b) `NEEDS EXTRACTION` and offers a fallback ("paired markers,
labelled as endpoints, do not interpolate") in case the per-step arrays cannot be recovered. This
module fixes which half is which, so that a later figure pass cannot quietly draw a smooth curve
through two points:

* the rank track HAS per-step values and they are asserted here;
* the collapse-evidence quantities have ONLY endpoints, and that absence is asserted too.

It also pins the statistic. `diag_d.py` computes `eff-rank` from its own inline formula, and that
formula is `RANK_VARIANTS["R3"]`, not `spectral.CANONICAL`. A test that only checked the numbers would
not have caught that; the two definitions are checked against each other on a matrix of known
spectrum below.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.spectral import RANK_VARIANTS, effective_rank
from morpheus.v2.research.rebase.p2 import p2_f8b_tracks


def test_vendored_logs_are_byte_identical_to_the_box_originals() -> None:
    """A log that has drifted is not the log the published numbers came from."""
    digests = p2_f8b_tracks.verify_sources()
    assert set(digests) == set(p2_f8b_tracks.SOURCES)
    for name, digest in digests.items():
        assert digest == p2_f8b_tracks.SOURCES[name][1]


def test_per_step_rank_track_is_recovered_at_the_recorded_steps() -> None:
    payload = p2_f8b_tracks.build()
    rows = payload["per_step_track"]["rows"]
    assert [row["step"] for row in rows] == [0, 25, 50, 100, 200, 400]
    # The published claim is "12.88 -> 1.00 by step 50", and it is a step-50 claim, not a step-400 one.
    assert rows[0]["rank_r3"] == pytest.approx(12.88)
    assert rows[2]["step"] == 50 and rows[2]["rank_r3"] == pytest.approx(1.00)
    assert all(row["rank_r3"] == pytest.approx(1.00) for row in rows[2:])
    # The collapse evidence co-measured on the same steps, which is what makes the track a diagnosis
    # rather than a rank curve: positives and worst negatives converge to within 1e-4.
    assert rows[2]["pos"] == pytest.approx(0.9993)
    assert rows[2]["worst_neg"] == pytest.approx(0.9993)
    assert rows[2]["min_margin"] == pytest.approx(-0.0001)


def test_the_track_statistic_is_R3_and_not_the_canonical_one() -> None:
    """`diag_d.py:50-51` is `(sum s)^2 / sum s^2` on centred, ROW-NORMALISED rows -- i.e. R3."""
    assert p2_f8b_tracks.RANK_STATISTIC == "R3"
    payload = p2_f8b_tracks.build()
    assert payload["rank_statistic"] == "R3"
    assert "not spectral.CANONICAL" in payload["rank_statistic_source"]

    # And the two are genuinely different functions, on a matrix whose spectrum is known by hand:
    # sigma proportional to (2, 1, 1) gives erank = 2*sqrt(2) under R1 and 8/3 under order 2.
    rng = np.random.default_rng(11)
    # Centre BEFORE orthonormalising, so the basis lies in the mean-zero subspace and the singular
    # values survive the centring `CANONICAL` applies. Orthonormalising first gives a basis with a
    # non-zero column mean, and the centred spectrum is then no longer (2, 1, 1).
    raw = rng.normal(size=(64, 3))
    basis = np.linalg.qr(raw - raw.mean(axis=0, keepdims=True))[0]
    x = basis * np.asarray([2.0, 1.0, 1.0])
    assert effective_rank(x, variant=RANK_VARIANTS["R1"]) == pytest.approx(2 * np.sqrt(2), rel=1e-9)
    assert effective_rank(x, variant=RANK_VARIANTS["R2"]) == pytest.approx(8 / 3, rel=1e-9)
    assert RANK_VARIANTS["R3"].order == 2 and RANK_VARIANTS["R3"].normalise_rows is True


def test_collapse_evidence_has_endpoints_only_and_carries_the_16_of_16_pinning() -> None:
    payload = p2_f8b_tracks.build()
    block = payload["endpoint_pairs_only"]
    assert block["per_step_array_retained"] is False
    assert "do not interpolate" in block["figure_instruction"]

    arms = block["arms"]
    assert set(arms) == {"A", "B", "C"}
    a = arms["A"]["endpoints"]
    # Every value F8(b) quotes for the collapsed arm, checked against the surviving log.
    assert a["WSI within-modality offdiag cos"] == {"before": 0.7089, "after": 0.9999}
    assert a["cross pos cos"] == {"before": 0.0538, "after": 0.9959}
    assert a["cross neg cos"] == {"before": 0.0816, "after": 0.9960}
    assert a["retrieval acc@1"]["before"] == 0.062 and a["retrieval acc@1"]["after"] == 0.000
    assert "chance 0.062" in a["retrieval acc@1"]["note"]
    # The panel's subject: a hard `matrix_rank` pinned at its structural maximum while everything
    # measured beside it collapses. Drawing this alone is what F8(b) exists to prevent.
    for arm in arms.values():
        pinned = arm["endpoints"]["z_biology matrix rank"]
        assert pinned["before"] == 16.0 and pinned["after"] == 16.0
        assert pinned["note"] == "max 16"


#: `v2/research/rebase/p2/figures/data/e0_run/` vendors the same two logs a second time, for the
#: figure-drawing pipeline. Two copies of one piece of evidence is the drift hazard this project
#: keeps re-learning, so if the second copy is present it must agree with the digest-enforced one.
DUPLICATES = {
    "collapse_diag.log": "figures/data/e0_run/collapse_diag.log",
    "diag_d.log": "figures/data/e0_run/d1_diag/diag_d.log",
}


def test_the_second_vendored_copy_of_each_log_agrees_with_this_one() -> None:
    """Compared on CONTENT, with line endings normalised.

    The `figures/data/` copies carry no `text` attribute, so a checkout with `core.autocrlf=true`
    rewrites their bytes and their sha256 changes while their content does not. Hashing them raw
    would turn every Windows checkout red for a non-problem; comparing normalised content still
    catches the failure that matters, which is the two copies genuinely diverging. That those copies
    are unprotected is a real gap and is recorded in
    `NOTEBOOK_ENTRIES/f8b_tracks_and_hill_order_inset_20260804T0630Z.md`, not fixed here -- the
    directory belongs to concurrent work.
    """
    root = p2_f8b_tracks.LOG_DIR.parent
    for name, relative in DUPLICATES.items():
        other = root / relative
        if not other.is_file():  # the figure pipeline is free not to vendor them
            continue
        mine = (p2_f8b_tracks.LOG_DIR / name).read_bytes().replace(b"\r\n", b"\n")
        theirs = other.read_bytes().replace(b"\r\n", b"\n")
        assert mine == theirs, f"{name} and {relative} have diverged"


def test_the_negative_pairs_are_the_higher_ones_in_the_collapsed_arm() -> None:
    """0.9960 > 0.9959: the reading that makes the instance evidence rather than an anomaly."""
    arms = p2_f8b_tracks.build()["endpoint_pairs_only"]["arms"]
    assert arms["A"]["endpoints"]["cross neg cos"]["after"] > arms["A"]["endpoints"]["cross pos cos"]["after"]
