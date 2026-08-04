"""Spike-recovery calibration — the CALIBRA methods contribution.

The problem this solves. Every confound-adjusted cross-modal result in this
field is reported without a sensitivity statement: when a paper residualises
site/cohort/purity out of an H&E->molecular analysis and finds a small effect,
nobody can say whether the adjustment *also* destroyed real signal, so a null is
uninterpretable. Three independent adversarial reviews of this project killed
three separate theses on exactly that objection.

The fix, borrowed in spirit from ERCC spike-ins in genomics and simulation-based
calibration in Bayesian methodology, and (to our knowledge) never applied to
certify a cross-modal morphology-molecular claim: inject a synthetic signal of
*known* strength into the molecular side, push it through the **identical**
pipeline including residualisation, and measure what comes out. That yields

  * a recovery curve  r_hat(r_true),
  * an empirical detection floor (smallest r_true recovered in >=80% of draws),
  * an attenuation slope,

against which the observed real-data effect can finally be read. It converts
"we adjusted and found nothing" into "our floor is 0.031 and we measured 0.068".

Spike construction. For a random unit direction ``v`` in molecular space and
``u`` in image space, let ``s = standardise(X u)`` and ``a = standardise(Y v)``.
We *replace* the v-component of Y with a signal having exactly the requested
correlation with s::

    a_perp = standardise(a - corr(s, a) * s)
    a_new  = r_true * s + sqrt(1 - r_true^2) * a_perp

so ``corr(s, a_new) == r_true`` by construction, before any residualisation.
Structured spikes (``direction="structured"``) use a real programme loading
vector for ``v`` instead of a random one, because a random direction is the
favourable case and would flatter the instrument.

TARGETED READOUT (2026-07-30 fix -- read this before changing the readout back).
The first implementation scored recovery with ``top_canonical_correlation``, a
*maximum* over ``n_components`` directions per side. The spike, however, lives on
exactly one known direction pair ``(u, v)``. On real data the ambient top-CCA sits
near 0.97, so a spike of r_true<=0.2 was invisible against it: every detection
floor came back NaN, level-0 recovery read 0.97 instead of ~0, and at r_true=0.2
the measured value *fell* (the replacement destroys pre-existing structure along
``v`` faster than the weak spike restores it). The instrument's headline output was
therefore unusable on the only data that matters.

The readout is now **direction-matched**: after pushing the spiked targets through
the identical residualisation, we score ``|corr(X_res u, Y_spiked_res v)|`` -- the
planted axis itself. The curve is then monotone in r_true and the floor is a real
detection threshold.

Level 0 is *not* zero on real data, and this is a genuine result rather than a
defect: ``a_new`` is built orthogonal to ``s``, but residualising two orthogonal
signals through a shared confound design induces a correlation between them
(~0.08-0.13 for the 99-column cancer+TSS design at n=2,530). How large that induced
floor is depends on how much of the drawn ``(u, v)`` lies in the design span, so it
varies draw to draw. The floor is therefore read with a **paired** test -- does a
draw beat its *own* level-0 value -- which is the point of reusing (u, v) across
levels within a draw.

Scale warning, CORRECTED 2026-08-04. The floor is in *single-direction* correlation
units. The **in-sample** ``observed`` (``top_canonical_correlation``) is a multivariate
maximum over ``n_components`` directions per side, is inflated by capacity -- its own
within-strata permutation null reaches p95 = 0.171 on the 99-column design -- and must
**never** be compared to the floor. That prohibition is real and unchanged.

What *was* wrong is the sentence that used to follow it. ``observed_matched_direction``
was offered as "the same-units comparator", and it is not a comparator at all: it is a
correlation along a **random** direction pair, which on a multi-column representation
carries essentially none of the channel and sits near zero however strong the channel
is, with a random sign on top. The flag built on it was a constant, not a measurement
(``NOTEBOOK_ENTRIES/observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md``).

The comparator that *is* on the floor's scale is a **held-out, fitted-direction**
statistic -- ``spectral.heldout_top_cca`` -- whose destroyed-pairing null is 0.06-0.09
against floors of 0.20-0.50, i.e. the capacity objection that condemns the in-sample
maximum is empirically false for the held-out one. The caller supplies it as
``SpikeRecoveryResult.channel_statistic`` and :func:`channel_clears_floor` grades it by
magnitude. See ``SUMMARY_SCHEMA_VERSION`` for what pre-fix artifacts carry, and
``direction_pairs`` / :func:`floors_from_recovery` for measuring the floor on the
channel's own fitted direction rather than a random one.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .residualise import cross_fitted_residuals
from .spectral import top_canonical_correlation

__all__ = ["SpikeRecoveryResult", "spike_targets", "spike_recovery_curve", "permutation_null",
           "floors_from_recovery", "channel_clears_floor", "SUMMARY_SCHEMA_VERSION"]

#: Schema version of :meth:`SpikeRecoveryResult.summary`'s dict, written into every
#: emitted row set so an artifact can say which contract it holds.
#:
#: * **1** — ``observed_above_floor`` compared the SIGNED ``observed_matched_direction``
#:   (a correlation along a RANDOM direction pair) against a positive
#:   ``detection_floor``. Both halves of that comparison were wrong, and the flag was
#:   in practice a **constant** rather than a measurement: for any multi-column
#:   representation a random direction pair carries ~none of the channel, so the
#:   left-hand side sits at ~0 whatever the channel does, and the sign made the
#:   verdict a coin flip on top of that. Every artifact emitted at schema 1 carries an
#:   ``observed_above_floor`` that must be ignored. Evidence:
#:   ``NOTEBOOK_ENTRIES/observed_above_floor_is_broken_and_every_channel_clears_20260804T2115Z.md``.
#: * **2** — the flag grades a caller-supplied **channel statistic** (a fitted-direction,
#:   out-of-fold quantity such as :func:`spectral.heldout_top_cca`) by **magnitude**
#:   against the floor, and is ``None`` with an explicit status when no channel
#:   statistic was supplied. It never invents a comparator from a random direction.
SUMMARY_SCHEMA_VERSION = 2

try:                                          # optional; the serial path is identical
    from joblib import Parallel, delayed
except ImportError:                           # pragma: no cover
    Parallel = None


def _map(fn, items, n_jobs: int):
    """Run ``fn(*item)`` over ``items``. Results are order-preserving and seed-fixed,
    so n_jobs changes wall-clock only -- never the numbers."""
    if n_jobs == 1 or Parallel is None or len(items) < 2:
        return [fn(*item) for item in items]
    return Parallel(n_jobs=n_jobs, prefer="processes")(delayed(fn)(*item) for item in items)


def permutation_null(x: np.ndarray, y: np.ndarray, design: np.ndarray, *, strata=None,
                     n_permutations: int = 100, n_components: int = 32, seed: int = 42,
                     n_jobs: int = 1) -> dict:
    """Top-CCA under destroyed cross-modal pairing — the chance level.

    Essential companion to the recovery curve. A top canonical correlation is a
    multivariate *maximum* over ``n_components`` directions per side, so at finite
    n it is substantially inflated by capacity alone: an observed 0.45 is
    meaningless until you know that chance gives, say, 0.42. Rows of ``y`` are
    permuted **within strata** (normally cancer type) so that cancer-level
    structure is preserved and only the patient-level pairing is destroyed —
    permuting globally would conflate the pairing with the cohort effect the
    residualisation is already meant to remove.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    rng = np.random.default_rng(seed)
    x_residual = cross_fitted_residuals(x, design, seed=seed)
    strata = np.zeros(len(x), dtype=int) if strata is None else np.asarray(strata)

    observed = top_canonical_correlation(x_residual, cross_fitted_residuals(y, design, seed=seed),
                                         n_components=n_components)
    # Draw every permutation order up front: identical numbers at any n_jobs.
    orders = []
    for _ in range(n_permutations):
        order = np.arange(len(y))
        for level in np.unique(strata):
            idx = np.flatnonzero(strata == level)
            order[idx] = rng.permutation(idx)
        orders.append((order,))

    def _one_permutation(order):
        permuted = cross_fitted_residuals(y[order], design, seed=seed)
        return top_canonical_correlation(x_residual, permuted, n_components=n_components)

    null = np.asarray(_map(_one_permutation, orders, n_jobs), dtype=np.float64)
    exceed = int(np.sum(null >= observed))
    return {
        "observed_top_cca": float(observed),
        "null_median": float(np.median(null)),
        "null_p95": float(np.percentile(null, 95)),
        "null_max": float(np.max(null)),
        "excess_over_null_median": float(observed - np.median(null)),
        "permutation_p": float((exceed + 1) / (n_permutations + 1)),
        "n_permutations": int(n_permutations),
    }


def _standardise(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=np.float64)
    centred = vector - vector.mean()
    scale = centred.std()
    return centred / scale if scale > 1e-12 else centred


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    """SIGNED Pearson r between two 1-D scores.

    Signed, not absolute, and that matters. The confound-induced level-0
    correlation has a random sign, so under |r| a draw whose induced baseline is
    negative gets *smaller* when the (positive) spike is added -- which destroys
    the paired comparison the floor depends on and reintroduces NaN floors. The
    spike is built with a known positive orientation (``a_new = r_true * s + ...``),
    so the signed statistic is the one that moves predictably with r_true.
    Magnitudes are taken at reporting time, never before the pairing.
    """
    a = _standardise(a)
    b = _standardise(b)
    if a.std() < 1e-12 or b.std() < 1e-12:
        return float("nan")
    return float(np.dot(a, b) / len(a))


def spike_targets(x: np.ndarray, y: np.ndarray, r_true: float, *, rng: np.random.Generator,
                  molecular_direction: np.ndarray | None = None,
                  image_direction: np.ndarray | None = None,
                  return_directions: bool = False):
    """Return a copy of ``y`` whose ``v``-component correlates with ``Xu`` at exactly ``r_true``.

    With ``return_directions=True`` returns ``(y_spiked, u, v)``. The caller needs
    ``(u, v)`` to score the spike on its own axis -- see the TARGETED READOUT note
    in the module docstring. The internal draw order (``u`` first, then ``v``) is
    part of the contract: tests reconstruct the directions by re-seeding an
    identical generator. **A supplied direction consumes no draw**, which is how
    ``molecular_direction`` has always behaved and is now also true of
    ``image_direction``.

    ``image_direction`` exists because until it did, no **direction-matched**
    floor could be measured. The floor is a statement about a direction pair, and
    with ``u`` drawn at random the floor was always measured on a pair that is
    random on at least the image side -- while the channel it grades
    (``spectral.heldout_top_cca``) uses a pair fitted on both sides. Pass
    ``spectral.heldout_cca_directions``' output here to plant on the axis the
    channel actually lives on. See
    ``NOTEBOOK_ENTRIES/PREDECLARED_direction_matched_floor_20260804T2230Z.md``
    for what that comparison does and does not settle: it matches the direction,
    it does **not** turn the oracle readout into an estimated one.
    """
    if not 0.0 <= r_true < 1.0:
        raise ValueError("r_true must lie in [0, 1)")
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if image_direction is None:
        u = rng.normal(size=x.shape[1])
    else:
        u = np.asarray(image_direction, dtype=np.float64).copy()
    norm_u = np.linalg.norm(u)
    if norm_u < 1e-12:
        raise ValueError("image_direction must be non-zero")
    u /= norm_u
    if molecular_direction is None:
        v = rng.normal(size=y.shape[1])
    else:
        v = np.asarray(molecular_direction, dtype=np.float64).copy()
    norm = np.linalg.norm(v)
    if norm < 1e-12:
        raise ValueError("molecular_direction must be non-zero")
    v /= norm

    s = _standardise(x @ u)
    a = _standardise(y @ v)
    rho = float(np.clip(np.dot(s, a) / len(s), -1.0, 1.0))
    a_perp = _standardise(a - rho * s)
    a_new = r_true * s + np.sqrt(max(0.0, 1.0 - r_true ** 2)) * a_perp
    # Replace, don't add: swap the v-component of Y for the constructed signal.
    # The update MUST be rescaled by the raw component's own s.d. ``a`` is
    # standardised while ``y @ v`` is not, so an unscaled ``a_new - a`` leaves a
    # residual (sigma - 1) * a in the target. That leftover carries the ambient
    # correlation rho, which then leaks into the readout: it put the level-0
    # baseline at 0.099 instead of ~0 and attenuated an r_true=0.6 spike to 0.27.
    # With the rescale, ``standardise(spiked @ v) == a_new`` exactly.
    sigma = float(np.std(y @ v))
    spiked = y + np.outer((a_new - a) * (sigma if sigma > 1e-12 else 1.0), v)
    return (spiked, u, v) if return_directions else spiked


def channel_clears_floor(channel_statistic: float, detection_floor: float) -> tuple:
    """Grade a channel statistic against a detection floor. Returns ``(verdict, status)``.

    ``verdict`` is ``True``/``False`` when the comparison can be made and ``None``
    when it cannot -- never a silent ``False``, because "we could not grade this"
    and "the channel is below its floor" are opposite scientific statements and
    the old flag conflated them for every artifact ever shipped.

    **The comparison is on MAGNITUDE.** A detection floor is a threshold on the
    size of an association; a channel is not less real for being negatively
    oriented, and no channel readout in this repository can even see the sign
    (``top_canonical_correlation``, ``heldout_top_cca`` and
    ``paired_absolute_correlation`` all take ``abs``). The signed convention of
    :func:`_correlation` is correct for the *paired* within-draw comparison it was
    written for and wrong here; carrying it into this unpaired comparison made
    three of four constructed sign-pairs return opposite verdicts for identical
    evidence.

    **``channel_statistic`` must come from the caller.** The floor grades *a
    channel*, and this module has no way to measure one: its own
    ``observed_matched_direction`` is a correlation along a **random** direction
    pair, which for any multi-column representation carries ~none of the channel
    and sits near zero whatever the channel does. Inventing the comparator here
    is what made the shipped flag a constant.
    """
    floor = float(detection_floor)
    channel = float(channel_statistic)
    if not np.isfinite(channel):
        return None, "ungraded_no_channel_statistic"
    if not np.isfinite(floor):
        return None, "ungraded_floor_did_not_resolve"
    return bool(abs(channel) > floor), "graded"


@dataclass
class SpikeRecoveryResult:
    levels: np.ndarray
    recovered: np.ndarray            # (n_levels, n_draws)
    detection_floor: float
    attenuation_slope: float
    observed: float = float("nan")
    n_components: int = 32
    recovery_fraction: float = 0.8
    meta: dict = field(default_factory=dict)

    delta: np.ndarray | None = None   # per-draw increment over that draw's level-0 baseline

    #: The channel this run's floor is to grade, supplied by the caller because
    #: this module cannot measure one (see :func:`channel_clears_floor`). It must
    #: be a FITTED-direction, out-of-fold statistic on the floor's single-direction
    #: scale -- ``spectral.heldout_top_cca`` is the one CALIBRA uses. Leave it NaN
    #: and the flag reports ``None`` with a status rather than a fabricated verdict.
    channel_statistic: float = float("nan")
    channel_statistic_name: str = ""

    def summary(self) -> dict:
        median = np.nanmedian(self.recovered, axis=1)
        delta = self.delta if self.delta is not None else self.recovered - self.recovered[0][None, :]
        verdict, status = channel_clears_floor(self.channel_statistic, self.detection_floor)
        return {
            "levels": self.levels.tolist(),
            "recovered_median": median.tolist(),
            "delta_median": np.nanmedian(delta, axis=1).tolist(),
            "delta_p10": np.nanpercentile(delta, 10, axis=1).tolist(),
            "recovered_p10": np.nanpercentile(self.recovered, 10, axis=1).tolist(),
            "recovered_p90": np.nanpercentile(self.recovered, 90, axis=1).tolist(),
            "detection_floor": self.detection_floor,
            "attenuation_slope": self.attenuation_slope,
            "observed": self.observed,
            # SCHEMA 2. At schema 1 this key compared the SIGNED
            # ``observed_matched_direction`` -- a random-direction correlation --
            # against a positive floor, and was a constant rather than a
            # measurement. It now grades the caller's channel statistic by
            # magnitude, and is None when there is nothing to grade.
            "summary_schema_version": SUMMARY_SCHEMA_VERSION,
            "observed_above_floor": verdict,
            "observed_above_floor_status": status,
            "channel_statistic": float(self.channel_statistic),
            "channel_statistic_name": self.channel_statistic_name,
            "floor_scale": "targeted_single_direction",
            "n_components": self.n_components,
            "recovery_fraction": self.recovery_fraction,
            **self.meta,
        }


def floors_from_recovery(levels: np.ndarray, recovered: np.ndarray, *,
                         recovery_fraction: float = 0.8) -> dict:
    """The two floors, read off a (n_levels, n_draws) recovery matrix.

    Lifted out of :func:`spike_recovery_curve` verbatim so that a caller scoring
    the spike with a **different readout** -- e.g. the channel's own
    ``spectral.heldout_top_cca``, which is what a like-for-like floor requires --
    reads its floor with the *same rule* rather than re-implementing it. The two
    floors are documented at length inside ``spike_recovery_curve``; the short
    version is that ``transmission_floor`` is paired (does a draw beat *its own*
    level 0) and answers the over-residualisation question, while
    ``detection_floor`` is unpaired (does it beat the level-0 p90) and is the
    conservative, quotable detection limit.
    """
    levels = np.asarray(levels, dtype=np.float64)
    recovered = np.asarray(recovered, dtype=np.float64)
    null_reference = (float(np.nanpercentile(recovered[0], 90))
                      if np.isfinite(recovered[0]).any() else 0.0)
    paired_hits = np.full(len(levels), np.nan)
    unpaired_hits = np.full(len(levels), np.nan)
    transmission_floor = float("nan")
    detection_floor = float("nan")
    for i, level in enumerate(levels):
        paired_hits[i] = np.nanmean(recovered[i] > recovered[0])
        unpaired_hits[i] = np.nanmean(recovered[i] > max(null_reference, 1e-9))
        if level <= 0:
            continue
        if np.isfinite(paired_hits[i]) and paired_hits[i] >= recovery_fraction \
                and not np.isfinite(transmission_floor):
            transmission_floor = float(level)
        if np.isfinite(unpaired_hits[i]) and unpaired_hits[i] >= recovery_fraction \
                and not np.isfinite(detection_floor):
            detection_floor = float(level)
    return {"detection_floor": detection_floor, "transmission_floor": transmission_floor,
            "null_reference_p90": null_reference,
            "paired_hit_rate": paired_hits.tolist(),
            "unpaired_hit_rate": unpaired_hits.tolist()}


def spike_recovery_curve(x: np.ndarray, y: np.ndarray, design: np.ndarray, *,
                         levels=(0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40),
                         n_draws: int = 25, n_components: int = 32, seed: int = 42,
                         recovery_fraction: float = 0.8,
                         molecular_directions: np.ndarray | None = None,
                         direction_pairs=None,
                         n_splits: int = 5, alpha: float = 1.0,
                         n_jobs: int = 1) -> SpikeRecoveryResult:
    """Push known-strength spikes through the full pipeline and measure recovery.

    ``design`` is the confound design matrix; residualisation is applied to the
    spiked data exactly as it is to the real data — that identity is the whole
    point, so do not "optimise" it away.

    ``direction_pairs`` pins BOTH sides of the planted axis: an iterable of
    ``(u, v)`` used one pair per draw, in order (cycled if shorter than
    ``n_draws``). With it the floor becomes **direction-matched** -- measured on
    the axis the channel actually occupies rather than on a random one, which is
    the comparison the shipped floors could not make. Mutually exclusive with
    ``molecular_directions``. Note what it does not fix: the readout still *knows*
    ``(u, v)`` while the channel's readout must estimate them, so a
    direction-matched oracle floor remains a lower bound on the floor that applies
    to a fitted-direction statistic. Use ``floors_from_recovery`` with a fitted
    readout for the like-for-like version.

    ``n_splits`` and ``alpha`` are the residualiser's own knobs, forwarded
    unchanged to :func:`cross_fitted_residuals`. They exist so that Track 2 can
    show the confound-induced level-0 baseline is a property of *shared
    residualisation* rather than of this particular fold count or shrinkage;
    the defaults reproduce every number produced before they were exposed.

    The detection floor is the smallest non-zero level whose recovered statistic
    exceeds the null (level 0) 90th percentile in at least ``recovery_fraction``
    of draws. Using the level-0 upper tail as the reference is what makes the
    floor a *detection* threshold rather than a bias estimate.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    levels = np.asarray(levels, dtype=np.float64)
    rng = np.random.default_rng(seed)

    resid = dict(n_splits=n_splits, alpha=alpha, seed=seed)
    x_residual = cross_fitted_residuals(x, design, **resid)
    y_residual = cross_fitted_residuals(y, design, **resid)

    if direction_pairs is not None and molecular_directions is not None:
        raise ValueError("pass direction_pairs OR molecular_directions, not both")

    # Draw every (direction, seed) up front so the result is identical regardless of
    # how the draws are scheduled across workers.
    plan = []
    for draw in range(n_draws):
        direction, image = None, None
        if direction_pairs is not None:
            if not len(direction_pairs):
                raise ValueError("direction_pairs must be non-empty")
            image, direction = direction_pairs[draw % len(direction_pairs)]
        elif molecular_directions is not None and len(molecular_directions):
            direction = molecular_directions[rng.integers(len(molecular_directions))]
        plan.append((int(rng.integers(1 << 31)), direction, image))

    def _one_draw(draw_seed: int, direction, image):
        """PAIRED: the same (u, v) is reused across every level within a draw."""
        column = np.full(len(levels), np.nan)
        matched = float("nan")
        for i, level in enumerate(levels):
            spiked, u, v = spike_targets(x, y, float(level),
                                         rng=np.random.default_rng(draw_seed),
                                         molecular_direction=direction,
                                         image_direction=image,
                                         return_directions=True)
            spiked_residual = cross_fitted_residuals(spiked, design, **resid)
            # TARGETED readout: score the planted axis, not a maximum over the whole
            # subspace. See the module docstring -- a max readout is swamped by
            # ambient structure and returns NaN floors on real data.
            column[i] = _correlation(x_residual @ u, spiked_residual @ v)
            if i == 0:
                # Same directions, same units, but on the UNSPIKED real targets.
                matched = _correlation(x_residual @ u, y_residual @ v)
        return column, matched

    results = _map(_one_draw, plan, n_jobs)
    recovered = np.column_stack([column for column, _ in results])
    observed_matched = float(np.nanmedian([m for _, m in results]))
    # The MAGNITUDE version. The signed median above is a median over draws whose
    # signs are random when (u, v) is random, so on a multi-column representation
    # it estimates approximately nothing; it is retained only because published
    # rows quote it. Neither is the channel -- see ``channel_clears_floor``.
    observed_matched_abs = float(np.nanmedian([abs(m) for _, m in results]))

    # PAIRED floor. The spike is built orthogonal to s, yet level 0 does NOT read
    # zero on real data: residualising two orthogonal signals through a shared
    # confound design induces a correlation between them (~0.08-0.13 for the 99-column
    # cancer+TSS design). That induced floor is draw-specific -- it depends on how much
    # of (u, v) lies in the design span -- so an unpaired threshold against the level-0
    # p90 is the wrong test: it charges every draw the WORST draw's induced baseline
    # and returns NaN even for spikes that are plainly recovered.
    #
    # Because the design is paired (same (u, v) at every level within a draw), the
    # correct test is a paired one: does this draw's recovered value exceed *its own*
    # level-0 value? The floor is the smallest level clearing that in
    # ``recovery_fraction`` of draws -- a paired sign test, which at 50% is chance and
    # at 80% is detection.
    #
    # TWO DIFFERENT FLOORS -- do not conflate them, and do not quote the paired one
    # as a detection limit.
    #
    #   transmission_floor (PAIRED): smallest r_true whose increment over the SAME
    #     draw's level-0 value is seen in >=recovery_fraction of draws. Pairing
    #     cancels the draw-specific induced baseline, so this is close to noiseless
    #     and typically resolves at the finest level on the grid. It answers "does
    #     the pipeline TRANSMIT a signal of this size, or destroy it?" -- the
    #     over-residualisation question. It is NOT a statement about statistical
    #     power, because a real analysis has no paired baseline to subtract.
    #
    #   detection_floor (UNPAIRED): smallest r_true whose recovered value clears the
    #     level-0 90th percentile in >=recovery_fraction of draws. This one carries
    #     the full draw-to-draw variability an actual single-shot analysis faces, so
    #     it is the conservative, quotable detection limit.
    #
    floors = floors_from_recovery(levels, recovered, recovery_fraction=recovery_fraction)
    null_reference = floors["null_reference_p90"]
    paired_hits = np.asarray(floors["paired_hit_rate"])
    unpaired_hits = np.asarray(floors["unpaired_hit_rate"])
    transmission_floor = floors["transmission_floor"]
    detection_floor = floors["detection_floor"]

    median_recovered = np.nanmedian(recovered, axis=1)
    finite = np.isfinite(levels) & np.isfinite(median_recovered)
    slope = float("nan")
    if finite.sum() >= 2:
        # d(recovered)/d(r_true): 1.0 means the adjustment cost nothing, 0 means it
        # destroyed the signal entirely. This is the attenuation the reviewers ask for.
        slope = float(np.polyfit(levels[finite], median_recovered[finite], 1)[0])

    delta = recovered - recovered[0][None, :]
    observed = top_canonical_correlation(x_residual, y_residual, n_components=n_components)

    # REAL-DATA GATE. With a targeted readout, level 0 must sit near the sampling
    # null (~1/sqrt(n)), not near 1. A baseline near 1 means the readout is picking
    # up ambient structure instead of the spike -- exactly the defect this fix
    # repairs -- and every floor below it would be meaningless.
    # Signed median sits near 0 (random sign); the MAGNITUDE is what the gate needs.
    baseline_signed = float(np.nanmedian(recovered[0]))
    baseline = float(np.nanmedian(np.abs(recovered[0])))
    null_scale = 3.0 / np.sqrt(max(len(x), 2))
    # The defect this gate exists to catch is a readout that reports AMBIENT structure
    # instead of the spike, whose signature is baseline ~= the multivariate top-CCA
    # (it was 0.97 against an ambient 0.89). Confound-induced baseline is a real but
    # much smaller effect and must not trip the gate. Hence a two-part bar: the
    # baseline must be either well under the ambient channel OR small in absolute
    # terms. The absolute allowance matters when the channel itself is weak, where
    # "half of ambient" would be an unreasonably tight bar for an effect that has
    # nothing to do with the defect. 0.15 is set just above what the 99-column
    # cancer+TSS design induces at n=2,530 (~0.08-0.13).
    gate_threshold = float(max(0.5 * float(observed), 0.15))
    return SpikeRecoveryResult(levels=levels, recovered=recovered, detection_floor=detection_floor,
                               attenuation_slope=slope, observed=float(observed),
                               n_components=n_components, recovery_fraction=recovery_fraction,
                               delta=delta,
                               meta={"null_reference_p90": null_reference, "n_draws": n_draws,
                                     "residualiser_n_splits": int(n_splits),
                                     "residualiser_alpha": float(alpha),
                                     "n_patients": int(len(x)),
                                     "baseline_recovered_median": baseline,
                                     "baseline_recovered_median_signed": baseline_signed,
                                     "baseline_null_scale": float(null_scale),
                                     "baseline_is_null_like": bool(baseline < gate_threshold),
                                     "baseline_gate_threshold": gate_threshold,
                                     "paired_hit_rate": paired_hits.tolist(),
                                     "unpaired_hit_rate": unpaired_hits.tolist(),
                                     "transmission_floor": transmission_floor,
                                     "confound_induced_baseline": baseline,
                                     "observed_matched_direction": observed_matched,
                                     "observed_matched_direction_abs": observed_matched_abs,
                                     "direction_source": ("fitted_pairs_supplied" if direction_pairs is not None
                                                          else "structured_molecular_direction"
                                                          if molecular_directions is not None else "random"),
                                     "observed_multivariate_top_cca": float(observed),
                                     "baseline_top_cca": float(top_canonical_correlation(
                                         x_residual, y_residual, n_components=n_components))})
