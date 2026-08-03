"""CALIBRA — Calibrated Cross-Modal Biological Recovery Audit.

The instrument that reports what effect size a confound-adjusted cross-modal
analysis *would have missed*. See `calibration.py` for the core contribution.
"""
from .calibration import SpikeRecoveryResult, spike_recovery_curve, spike_targets
from .residualise import confound_design, cross_fitted_residuals
from .spectral import (CANONICAL, RANK_VARIANTS, RankVariant, cca_spectrum, effective_rank,
                       effective_rank_all_variants, top_canonical_correlation)

__all__ = [
    "SpikeRecoveryResult", "spike_recovery_curve", "spike_targets",
    "confound_design", "cross_fitted_residuals",
    "cca_spectrum", "effective_rank", "effective_rank_all_variants",
    "top_canonical_correlation", "CANONICAL", "RANK_VARIANTS", "RankVariant",
]
