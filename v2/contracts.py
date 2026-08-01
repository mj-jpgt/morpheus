"""Versioned contracts for V2 checkpoints and frozen representations.

The first V2 run exported heads that had never occurred in a loss.  This
module makes that state impossible to silently evaluate again: an artifact
must name the states it trained, and consumers must respect that declaration.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


# Version 4 adds source provenance and V2.1 teacher-mode metadata.  Readers
# intentionally accept v3 baseline artifacts so the complete historical suite
# remains comparable under one evaluation protocol.
ARTIFACT_VERSION = 4
IDENTITY_STATES = frozenset(("wsi_identity", "rna_identity", "full_identity"))
BIOLOGY_STATES = frozenset(("wsi_biology", "rna_biology", "full_biology"))
LEGACY_IDENTITY_ALIASES = {"wsi_identity": "wsi", "rna_identity": "rna"}


def declared_states(raw: np.lib.npyio.NpzFile) -> frozenset[str]:
    """Return explicitly trained states; legacy files only expose identities.

    Legacy exports have no trustworthy declaration.  Their identity views were
    directly contrastively trained, while their patient and biology states were
    not consistently supervised, so they are intentionally withheld.
    """
    if "trained_states" in raw.files:
        return frozenset(str(item) for item in raw["trained_states"].astype(str))
    states = {key for key in IDENTITY_STATES if key in raw.files}
    states.update(state for state, legacy in LEGACY_IDENTITY_ALIASES.items() if legacy in raw.files)
    return frozenset(states)


def require_state(raw: np.lib.npyio.NpzFile, state: str) -> np.ndarray | None:
    """Return a declared state or ``None``; never guess a fallback view."""
    if state not in declared_states(raw):
        return None
    key = state if state in raw.files else LEGACY_IDENTITY_ALIASES.get(state)
    return raw[key].astype(np.float32) if key is not None and key in raw.files else None


def validate_artifact(path: str | Path) -> dict[str, object]:
    """Validate the minimal frozen-representation schema without mutation."""
    with np.load(path, allow_pickle=False) as raw:
        required = {"patient_ids", "cancers", "split"}
        missing = required - set(raw.files)
        if missing:
            raise ValueError(f"artifact is missing {sorted(missing)}")
        rows = len(raw["patient_ids"])
        if any(len(raw[key]) != rows for key in required):
            raise ValueError("artifact identifier arrays have inconsistent lengths")
        states = declared_states(raw)
        if "artifact_version" in raw.files and int(raw["artifact_version"].item()) > ARTIFACT_VERSION:
            raise ValueError("artifact was produced by a newer unsupported contract")
        missing_states = {state for state in states if state not in raw.files and LEGACY_IDENTITY_ALIASES.get(state) not in raw.files}
        if missing_states:
            raise ValueError(f"artifact declares missing states {sorted(missing_states)}")
        for state in states:
            key = state if state in raw.files else LEGACY_IDENTITY_ALIASES[state]
            value = raw[key]
            if value.ndim != 2 or value.shape[0] != rows or not np.isfinite(value).all():
                raise ValueError(f"invalid representation state {state!r}")
    return {"artifact_version": ARTIFACT_VERSION, "n_patients": rows, "trained_states": sorted(states),
            "has_manifest": "manifest_json" in raw.files}


def trainable_state_names(*names: str | None) -> tuple[str, ...]:
    """Stable, duplicate-free state declaration for manifests and NPZ files."""
    return tuple(dict.fromkeys(name for name in names if name))


def assert_module_used(module, loss) -> None:
    """Fail fast if a state head did not participate in the supplied loss.

    This is called by a small unit test and deliberately examines gradients
    after ``backward``.  It is not part of the performance-critical loop.
    """
    del loss
    if not any(parameter.requires_grad and parameter.grad is not None and parameter.grad.detach().abs().sum() > 0
               for parameter in module.parameters()):
        raise AssertionError("exported state head received no gradient from its declared objective")
