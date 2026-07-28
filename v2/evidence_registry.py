"""Closed, versioned evidence registry for MORPHEUS hypothesis cards.

The registry is deliberately local and data-free: it contains reviewed
mechanism/assay facts, never patient labels, paired RNA, outcomes, or web
search results.  Qwen may cite only chunks selected from this registry.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .hypothesis_cards import ASSAYS, MECHANISM_BY_NAME, PROGRAMMES


REGISTRY_VERSION = "morpheus_closed_evidence_registry_v1"


@dataclass(frozen=True)
class EvidenceChunk:
    """A reviewed, local evidence unit that may be cited by a generated card."""

    id: str
    title: str
    text: str
    supported_programmes: tuple[str, ...]
    permitted_assays: tuple[str, ...]
    source_kind: str = "curated_ontology"

    def validate(self) -> None:
        if not self.id or not self.title or not self.text:
            raise ValueError("evidence chunks need non-empty id, title, and text")
        if set(self.supported_programmes) - PROGRAMMES:
            raise ValueError(f"evidence {self.id!r} has unsupported programmes")
        if set(self.permitted_assays) - ASSAYS:
            raise ValueError(f"evidence {self.id!r} has unsupported assays")
        if self.source_kind not in {"curated_ontology", "reviewed_local_literature"}:
            raise ValueError("unsupported evidence source_kind")


class EvidenceRegistry:
    """Immutable index that exposes only local, pre-reviewed evidence chunks."""

    def __init__(self, chunks: Iterable[EvidenceChunk], *, version: str = REGISTRY_VERSION):
        self.version = str(version)
        values = tuple(chunks)
        if not values:
            raise ValueError("evidence registry cannot be empty")
        if len({chunk.id for chunk in values}) != len(values):
            raise ValueError("evidence registry IDs must be unique")
        for chunk in values:
            chunk.validate()
        self._chunks = {chunk.id: chunk for chunk in values}
        canonical = {"version": self.version, "chunks": [asdict(c) for c in sorted(values, key=lambda c: c.id)]}
        self.digest = hashlib.sha256(json.dumps(canonical, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    @property
    def chunks(self) -> tuple[EvidenceChunk, ...]:
        return tuple(self._chunks[key] for key in sorted(self._chunks))

    def get(self, evidence_id: str) -> EvidenceChunk:
        try:
            return self._chunks[str(evidence_id)]
        except KeyError as exc:
            raise KeyError(f"unknown local evidence ID {evidence_id!r}") from exc

    def retrieve(
        self,
        programmes: Sequence[str],
        *,
        assay: str | None = None,
        limit: int = 4,
    ) -> tuple[EvidenceChunk, ...]:
        """Deterministically retrieve local evidence by declared programme overlap.

        This intentionally has no network or embedding-model dependency.  It is
        a transparent closed-RAG baseline, and cannot leak test-patient labels.
        """
        if limit < 1:
            raise ValueError("limit must be positive")
        requested = set(map(str, programmes)) & PROGRAMMES
        if assay is not None and assay not in ASSAYS:
            raise ValueError("unsupported assay")
        scored: list[tuple[int, str, EvidenceChunk]] = []
        for chunk in self._chunks.values():
            overlap = len(requested & set(chunk.supported_programmes))
            assay_bonus = 1 if assay is not None and assay in chunk.permitted_assays else 0
            if overlap or assay_bonus:
                scored.append((overlap * 2 + assay_bonus, chunk.id, chunk))
        return tuple(item[2] for item in sorted(scored, key=lambda x: (-x[0], x[1]))[:limit])

    def to_json(self, path: str | Path) -> None:
        payload = {"version": self.version, "digest": self.digest, "chunks": [asdict(c) for c in self.chunks]}
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    @classmethod
    def from_json(cls, path: str | Path) -> "EvidenceRegistry":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping) or set(payload) - {"version", "digest", "chunks"}:
            raise ValueError("invalid evidence registry JSON")
        chunks = []
        for raw in payload.get("chunks", []):
            if set(raw) != {"id", "title", "text", "supported_programmes", "permitted_assays", "source_kind"}:
                raise ValueError("invalid evidence chunk JSON")
            chunks.append(EvidenceChunk(
                id=str(raw["id"]), title=str(raw["title"]), text=str(raw["text"]),
                supported_programmes=tuple(map(str, raw["supported_programmes"])),
                permitted_assays=tuple(map(str, raw["permitted_assays"])), source_kind=str(raw["source_kind"]),
            ))
        registry = cls(chunks, version=str(payload.get("version", REGISTRY_VERSION)))
        if "digest" in payload and str(payload["digest"]) != registry.digest:
            raise ValueError("evidence registry digest does not match its content")
        return registry


def default_evidence_registry() -> EvidenceRegistry:
    """Return the finite reviewed ontology shipped with MORPHEUS.

    These are assay-falsification templates, not clinical guidance or evidence
    about any individual patient.
    """
    chunks = []
    for name, mechanism in sorted(MECHANISM_BY_NAME.items()):
        programmes = tuple(programme for programme, _ in mechanism.required)
        chunks.append(EvidenceChunk(
            id=f"mechanism:{name}", title=name.replace("_", " "),
            text=("Predeclared research mechanism template. " + mechanism.predicted_observation),
            supported_programmes=programmes, permitted_assays=(mechanism.assay,),
        ))
    return EvidenceRegistry(chunks)
