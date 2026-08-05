"""Guards for the multi-partition validation of the inductive channel result.

Predeclared in ``NOTEBOOK_ENTRIES/PREDECLARED_inductive_channel_split_stability_20260805T0015Z.md``.
That run re-measures ``d2_h::wsi_biology``'s channel retention on eight discovery/exposure
partitions instead of the one the 23:45 entry reported 0.9966 on, and its whole value depends on
two premises that are checked here rather than assumed:

* **the partitions are actually different.** ``--split-seed`` is threaded only into
  ``exposure_split``; if it were not wired, or if ``exposure_split`` ignored its seed, eight runs
  would return eight copies of one number and a spread of 0.000 would be reported as stability.
  ``test_two_split_seeds_produce_genuinely_different_partitions`` is §5.1 of the predeclaration:
  same fold sizes, every cancer still on both sides, and an exposure-set overlap near the 1/3 two
  independent halves of a cohort must have -- never near 1.
* **the coverage number the run regresses retention against is the operator's own.** The driver
  reports how many exposure rows carry a site the discovery fold saw ``min_site_count`` times, and
  that count must come from the same ``SitePooling`` the operator adjusts with, not from a second
  count computed beside it.

``test_split_seed_defaults_to_the_published_command`` pins the default-preserving property: the
published command carries no ``--split-seed``, so the flag's default must resolve to ``--seed``.
"""
from __future__ import annotations

import numpy as np
import pytest

from morpheus.v2.calibra.inductive_adjustment import ConfoundAdjustmentOperator
from morpheus.v2.research.rebase.nature.p1_evidence import inductive_channel
from morpheus.v2.research.rebase.nature.p4_certification.p4_certify import exposure_split

#: The partitions the predeclaration commits to, in its order.
PREDECLARED_SPLIT_SEEDS = (42, 7, 11, 23, 101, 555, 2718, 31337)

SITES = ["01", "02", "05", "06", "09", "13", "17", "22", "31", "44"]
CANCERS = ["BRCA", "LUAD", "COAD", "GBM", "STAD", "KIRC"]


def _cohort(n: int = 600, seed: int = 0):
    rng = np.random.default_rng(seed)
    site = rng.choice(SITES, size=n)
    cancers = rng.choice(CANCERS, size=n)
    patient_ids = np.asarray([f"TCGA-{s}-{i:04d}" for i, s in enumerate(site)])
    return patient_ids, cancers


def _exposure_ids(patient_ids, cancers, seed, fraction=0.5):
    discovery = exposure_split(cancers, discovery_fraction=fraction, seed=seed)
    return set(patient_ids[~discovery].tolist()), set(patient_ids[discovery].tolist())


# --- §5.1: the partitions are actually different ------------------------------------------

@pytest.mark.parametrize("other", PREDECLARED_SPLIT_SEEDS[1:])
def test_two_split_seeds_produce_genuinely_different_partitions(other):
    patient_ids, cancers = _cohort()
    baseline, _ = _exposure_ids(patient_ids, cancers, 42)
    alternative, _ = _exposure_ids(patient_ids, cancers, other)
    assert len(alternative) == len(baseline)            # matched n, or the arms are not comparable
    assert alternative != baseline
    jaccard = len(baseline & alternative) / len(baseline | alternative)
    # Two independent halves of the same cohort share ~1/3 of their union. The predeclaration
    # voids the run above 0.9; the band here is wide enough not to be a flake and narrow enough
    # to catch a seed that is not wired through.
    assert 0.15 < jaccard < 0.55, jaccard


def test_every_predeclared_split_seed_keeps_every_cancer_on_both_sides():
    """``on_unseen_level='refuse'`` stops a run whose exposure fold carries an unseen cancer."""
    patient_ids, cancers = _cohort()
    for seed in PREDECLARED_SPLIT_SEEDS:
        exposure, discovery = _exposure_ids(patient_ids, cancers, seed)
        exposure_cancers = {c for i, c in zip(patient_ids, cancers) if i in exposure}
        discovery_cancers = {c for i, c in zip(patient_ids, cancers) if i in discovery}
        assert exposure_cancers == discovery_cancers == set(CANCERS), seed


def test_the_eight_partitions_are_eight_distinct_partitions():
    patient_ids, cancers = _cohort()
    digests = {frozenset(_exposure_ids(patient_ids, cancers, seed)[0])
               for seed in PREDECLARED_SPLIT_SEEDS}
    assert len(digests) == len(PREDECLARED_SPLIT_SEEDS)


# --- the coverage covariate is the operator's own count -----------------------------------

def test_site_coverage_is_read_from_the_operators_own_pooling_rule():
    """The count the run regresses retention against must be the operator's, not a copy of it."""
    patient_ids, cancers = _cohort(n=400, seed=3)
    discovery = exposure_split(cancers, discovery_fraction=0.5, seed=42)
    fit_rows, score_rows = np.flatnonzero(discovery), np.flatnonzero(~discovery)
    rng = np.random.default_rng(0)
    matrix = rng.normal(size=(len(patient_ids), 6))
    operator = ConfoundAdjustmentOperator.fit(
        matrix[fit_rows], {"cancer": cancers[fit_rows]}, ["cancer", "tss"],
        patient_ids=patient_ids[fit_rows], site_column="tss", min_site_count=10)

    pooled, report = operator.site_pooling.apply(patient_ids[score_rows])
    covered = int(len(score_rows) - report["n_pooled_to_other"])
    assert covered == int(np.sum(pooled != "OTHER"))
    # and the rule really is the *discovery* fold's, evaluated at a count of zero for a site the
    # discovery fold never saw
    assert set(pooled.tolist()) - {"OTHER"} <= set(operator.site_pooling.frequent)


def test_split_seed_defaults_to_the_published_command():
    """``--split-seed`` absent must mean ``--seed``, or every published number moves."""
    resolve = inductive_channel.resolve_split_seed
    assert resolve(42, -1) == 42          # the published command, unmoved
    assert resolve(7, -1) == 7
    assert resolve(42, 7) == 7            # and set, it is the split seed and only the split seed
    assert resolve(42, 31337) == 31337


def test_the_default_split_seed_reproduces_the_published_partition_exactly():
    patient_ids, cancers = _cohort()
    published, _ = _exposure_ids(patient_ids, cancers, 42)
    default, _ = _exposure_ids(patient_ids, cancers,
                               inductive_channel.resolve_split_seed(42, -1))
    assert default == published
