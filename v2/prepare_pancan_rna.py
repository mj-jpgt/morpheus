"""Prepare the wide TCGA PanCan EBPlusPlus RNA matrix for MORPHEUS.

The TCGA PanCan RNA release is a gene-by-sample TSV (rather than the
patient-by-gene table used by :mod:`discovery_targets`).  Loading and
transposing the roughly 2 GB source in one pandas call is unnecessarily
fragile on a shared machine.  This module makes two streaming passes over the
source and uses disk-backed arrays while reducing it:

* source sample barcodes are mapped to canonical 12-character participants;
* all source samples belonging to a participant are averaged deterministically;
* ``SYMBOL|ENTREZ``/``ENSG...|SYMBOL`` labels are normalised to symbols;
* repeated gene symbols are averaged using all observed source values;
* a patient-by-gene float32 parquet table and a complete provenance manifest
  are written atomically.

It is deliberately strict about malformed rows, duplicate *sample columns*,
and genuinely ambiguous gene labels.  Unmapped labels (for example ``?|0``)
are excluded because they cannot be scored against a symbol-based GMT; the
number and examples are retained in the manifest.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import Iterable

import numpy as np
import pandas as pd

from .canonical_registry import RegistryError, canonical_tcga_patient_id


class RNAPreparationError(ValueError):
    """Raised when a source matrix cannot be converted without guessing."""


_GENE_SYMBOL = re.compile(r"^[A-Z][A-Z0-9_.-]*$")
_ENSEMBL = re.compile(r"^ENS(?:G|T)[0-9]+(?:\.[0-9]+)?$", re.IGNORECASE)
_UNMAPPED = {"", "?", "NA", "N/A", "NAN", "NULL", "NONE", "-"}


@dataclass(frozen=True)
class SourceLayout:
    sample_ids: tuple[str, ...]
    patient_ids: tuple[str, ...]
    patient_for_sample: tuple[str, ...]
    genes: tuple[str, ...]
    source_gene_rows: int
    dropped_gene_rows: int
    dropped_gene_examples: tuple[str, ...]
    source_sha256: str


def _clean_tsv_scalar(value: str) -> str:
    """Remove the simple RFC4180 quoting used by the PanCan release header."""
    text = str(value).strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        # A quoted TSV field represents literal quotes as doubled quotes.
        return text[1:-1].replace('""', '"').strip()
    return text


def _normalise_gene_label(value: str, *, line_number: int) -> str | None:
    """Return one unambiguous HGNC-style symbol, or ``None`` when unmapped.

    TCGA normally uses ``A1BG|1``.  Some mirrors expose
    ``ENSG...|A1BG|1``.  If two plausible symbols are present, choosing either
    would silently alter pathway scores, so that case is a hard error.
    """
    text = _clean_tsv_scalar(value).upper()
    if text in _UNMAPPED:
        return None
    pieces = [piece.strip().upper() for piece in text.split("|")]
    candidates = [
        piece for piece in pieces
        if piece not in _UNMAPPED
        and not piece.isdigit()
        and not _ENSEMBL.fullmatch(piece)
        and _GENE_SYMBOL.fullmatch(piece) is not None
    ]
    candidates = list(dict.fromkeys(candidates))
    if not candidates:
        return None
    if len(candidates) != 1:
        raise RNAPreparationError(
            f"ambiguous gene label at source line {line_number}: {value!r} "
            f"maps to {candidates}"
        )
    return candidates[0]


def _read_header_and_layout(source: Path) -> SourceLayout:
    """First streaming pass: validate TSV shape, canonicalise identifiers, hash."""
    digest = sha256()
    source_gene_rows = dropped = 0
    dropped_examples: list[str] = []
    gene_counter: Counter[str] = Counter()
    with source.open("rb") as handle:
        first = handle.readline()
        if not first:
            raise RNAPreparationError(f"RNA source is empty: {source}")
        digest.update(first)
        try:
            header = first.decode("utf-8-sig").rstrip("\r\n").split("\t")
        except UnicodeDecodeError as error:
            raise RNAPreparationError(f"RNA source header is not UTF-8: {source}") from error
        if len(header) < 2 or not header[0].strip():
            raise RNAPreparationError("RNA TSV needs a gene-ID column and at least one sample column")
        sample_ids = tuple(_clean_tsv_scalar(item) for item in header[1:])
        if any(not item for item in sample_ids):
            raise RNAPreparationError("RNA TSV has an empty sample column name")
        duplicate_samples = [name for name, count in Counter(sample_ids).items() if count > 1]
        if duplicate_samples:
            raise RNAPreparationError(
                "RNA TSV has duplicate source sample columns; cannot safely distinguish them: "
                f"{duplicate_samples[:5]}"
            )
        try:
            patient_for_sample = tuple(canonical_tcga_patient_id(sample) for sample in sample_ids)
        except RegistryError as error:
            raise RNAPreparationError(
                "Every RNA matrix column must be a TCGA patient/sample barcode; "
                f"received an invalid header value. {error}"
            ) from error
        patient_ids = tuple(sorted(set(patient_for_sample)))
        expected_tabs = len(sample_ids)
        for line_number, raw in enumerate(handle, start=2):
            digest.update(raw)
            # This release is a plain numeric TSV.  Rejecting irregular rows
            # here avoids pandas' permissive column-shift behaviour later.
            if raw.count(b"\t") != expected_tabs:
                raise RNAPreparationError(
                    f"malformed source row {line_number}: expected {expected_tabs + 1} fields"
                )
            try:
                label = _clean_tsv_scalar(raw.split(b"\t", 1)[0].decode("utf-8"))
            except UnicodeDecodeError as error:
                raise RNAPreparationError(f"gene label at source line {line_number} is not UTF-8") from error
            symbol = _normalise_gene_label(label, line_number=line_number)
            source_gene_rows += 1
            if symbol is None:
                dropped += 1
                if len(dropped_examples) < 10:
                    dropped_examples.append(label)
            else:
                gene_counter[symbol] += 1
    if not gene_counter:
        raise RNAPreparationError("no mappable gene symbols were found in the RNA source")
    return SourceLayout(
        sample_ids=sample_ids,
        patient_ids=patient_ids,
        patient_for_sample=patient_for_sample,
        genes=tuple(sorted(gene_counter)),
        source_gene_rows=source_gene_rows,
        dropped_gene_rows=dropped,
        dropped_gene_examples=tuple(dropped_examples),
        source_sha256=digest.hexdigest(),
    )


def _sample_grouping(layout: SourceLayout) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return sorting/reduction arrays to aggregate source samples per patient."""
    patient_index = {patient: index for index, patient in enumerate(layout.patient_ids)}
    membership = np.asarray([patient_index[patient] for patient in layout.patient_for_sample], dtype=np.int64)
    order = np.argsort(membership, kind="stable")
    counts = np.bincount(membership, minlength=len(layout.patient_ids))
    starts = np.r_[0, np.cumsum(counts)[:-1]].astype(np.int64, copy=False)
    return order, starts, counts


def _aggregate_gene_values(values: np.ndarray, order: np.ndarray, starts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sum finite values and counts within each canonical patient."""
    ordered = values[order]
    finite = np.isfinite(ordered)
    sums = np.add.reduceat(np.where(finite, ordered, 0.0), starts).astype(np.float32, copy=False)
    counts = np.add.reduceat(finite.astype(np.uint32), starts).astype(np.uint16, copy=False)
    return sums, counts


def _iter_source_chunks(source: Path, *, chunk_genes: int, sample_ids: tuple[str, ...]) -> Iterable[pd.DataFrame]:
    try:
        yield from pd.read_csv(
            source,
            sep="\t",
            index_col=0,
            chunksize=chunk_genes,
            # Passing a scalar dtype also targets the index column in some
            # pandas/C-parser versions.  A per-sample mapping keeps the gene
            # identifier as text while retaining C-level float parsing.
            dtype={sample: np.float32 for sample in sample_ids},
            na_values=["", "NA", "N/A", "NaN", "nan"],
            keep_default_na=True,
            low_memory=False,
        )
    except (ValueError, pd.errors.ParserError) as error:
        raise RNAPreparationError(f"could not parse numeric RNA values from {source}: {error}") from error


def _sidecar_path(output: Path, suffix: str) -> Path:
    return output.with_name(f"{output.stem}{suffix}")


def prepare_pancan_rna(
    source: str | Path,
    output: str | Path,
    *,
    chunk_genes: int = 128,
    patient_chunk: int = 512,
    work_dir: str | Path | None = None,
    overwrite: bool = False,
    keep_work: bool = False,
) -> dict[str, object]:
    """Convert a wide TCGA RNA matrix to patient-by-gene parquet.

    ``chunk_genes`` controls source parsing RAM.  The dense intermediate is a
    memmap (roughly 6 bytes per patient-gene cell: float32 sum + uint16 count),
    not host RAM.  ``output`` has a ``.manifest.json`` and ``.sample_map.parquet``
    sidecar suitable for the discovery-run provenance audit.
    """
    source_path, output_path = Path(source), Path(output)
    if chunk_genes < 1 or patient_chunk < 1:
        raise RNAPreparationError("chunk sizes must be positive")
    if not source_path.is_file():
        raise RNAPreparationError(f"RNA source does not exist: {source_path}")
    if output_path.exists() and not overwrite:
        raise RNAPreparationError(f"output already exists (use --overwrite): {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    layout = _read_header_and_layout(source_path)
    gene_index = {gene: index for index, gene in enumerate(layout.genes)}
    order, starts, _ = _sample_grouping(layout)
    temporary_root = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="morpheus-rna-", dir=output_path.parent))
    created_temp = work_dir is None
    temporary_root.mkdir(parents=True, exist_ok=True)
    sums_path, counts_path = temporary_root / "gene_sums.f32.mmap", temporary_root / "gene_counts.u16.mmap"
    parquet_tmp = output_path.with_suffix(output_path.suffix + ".tmp")
    sample_map_path = _sidecar_path(output_path, ".sample_map.parquet")
    manifest_path = _sidecar_path(output_path, ".manifest.json")
    if parquet_tmp.exists():
        parquet_tmp.unlink()
    try:
        shape = (len(layout.patient_ids), len(layout.genes))
        sums = np.memmap(sums_path, dtype=np.float32, mode="w+", shape=shape)
        counts = np.memmap(counts_path, dtype=np.uint16, mode="w+", shape=shape)
        sums[:] = 0.0; counts[:] = 0
        duplicate_rows: Counter[str] = Counter()
        source_rows_seen = 0
        for chunk_offset, chunk in enumerate(
            _iter_source_chunks(source_path, chunk_genes=chunk_genes, sample_ids=layout.sample_ids)
        ):
            if tuple(chunk.columns.astype(str)) != layout.sample_ids:
                raise RNAPreparationError("RNA parser changed source sample columns after preflight")
            values = chunk.to_numpy(dtype=np.float32, copy=False)
            if values.shape[1] != len(layout.sample_ids):
                raise RNAPreparationError("RNA parser returned an unexpected number of sample columns")
            for local_row, raw_label in enumerate(chunk.index.astype(str)):
                line_number = 2 + chunk_offset * chunk_genes + local_row
                symbol = _normalise_gene_label(raw_label, line_number=line_number)
                source_rows_seen += 1
                if symbol is None:
                    continue
                gene = gene_index.get(symbol)
                if gene is None:
                    raise RNAPreparationError(f"gene {symbol} was absent from first-pass index")
                row_sum, row_count = _aggregate_gene_values(values[local_row], order, starts)
                # uint16 overflow is a data integrity failure, not a reason to
                # wrap counts and corrupt duplicate-gene reduction.
                if np.any(counts[:, gene].astype(np.uint32) + row_count.astype(np.uint32) > np.iinfo(np.uint16).max):
                    raise RNAPreparationError(f"too many observed values while collapsing duplicate gene {symbol}")
                sums[:, gene] += row_sum
                counts[:, gene] += row_count
                duplicate_rows[symbol] += 1
        if source_rows_seen != layout.source_gene_rows:
            raise RNAPreparationError(
                f"source row count changed between passes ({layout.source_gene_rows} -> {source_rows_seen})"
            )
        sums.flush(); counts.flush()
        import pyarrow as pa
        import pyarrow.parquet as pq
        writer: pq.ParquetWriter | None = None
        try:
            for start in range(0, shape[0], patient_chunk):
                stop = min(start + patient_chunk, shape[0])
                local_counts = np.asarray(counts[start:stop], dtype=np.float32)
                local_sums = np.asarray(sums[start:stop], dtype=np.float32)
                matrix = np.divide(local_sums, local_counts, out=np.full_like(local_sums, np.nan), where=local_counts > 0)
                frame = pd.DataFrame(matrix, columns=layout.genes, copy=False)
                frame.insert(0, "patient_id", list(layout.patient_ids[start:stop]))
                table = pa.Table.from_pandas(frame, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(parquet_tmp, table.schema, compression="zstd")
                writer.write_table(table)
        finally:
            if writer is not None:
                writer.close()
        if writer is None:
            raise RNAPreparationError("no patient rows were written")
        os.replace(parquet_tmp, output_path)
        sample_map = pd.DataFrame({"source_sample_id": layout.sample_ids, "patient_id": layout.patient_for_sample})
        sample_map.to_parquet(sample_map_path, index=False)
        repeated_genes = {gene: count for gene, count in duplicate_rows.items() if count > 1}
        manifest: dict[str, object] = {
            "schema_version": "1.0",
            "source": str(source_path.resolve()),
            "source_sha256": layout.source_sha256,
            "source_gene_rows": layout.source_gene_rows,
            "mappable_gene_rows": layout.source_gene_rows - layout.dropped_gene_rows,
            "dropped_unmapped_gene_rows": layout.dropped_gene_rows,
            "dropped_unmapped_gene_examples": list(layout.dropped_gene_examples),
            "unique_gene_symbols": len(layout.genes),
            "duplicate_gene_symbols_collapsed": len(repeated_genes),
            "duplicate_gene_row_count": int(sum(count - 1 for count in repeated_genes.values())),
            "duplicate_gene_examples": dict(list(sorted(repeated_genes.items()))[:20]),
            "source_sample_count": len(layout.sample_ids),
            "canonical_patient_count": len(layout.patient_ids),
            "multi_sample_patient_count": int(sum(count > 1 for count in Counter(layout.patient_for_sample).values())),
            "max_samples_per_patient": int(max(Counter(layout.patient_for_sample).values())),
            "sample_header_sha256": sha256("\n".join(layout.sample_ids).encode("utf-8")).hexdigest(),
            "patient_map_sha256": sha256(sample_map.to_csv(index=False, lineterminator="\n").encode("utf-8")).hexdigest(),
            "output": str(output_path.resolve()),
            "sample_map": str(sample_map_path.resolve()),
            "dtype": "float32",
            "reduction": "mean over all finite source values for each canonical patient and gene symbol",
            "chunk_genes": chunk_genes,
            "patient_chunk": patient_chunk,
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest
    finally:
        # Explicitly release Windows file handles before cleanup.
        for name in ("sums", "counts"):
            value = locals().get(name)
            if isinstance(value, np.memmap):
                value.flush()
                # numpy does not provide a public close() for memmap, but the
                # backing mmap must be closed before a Windows scratch file
                # can be removed.  Guard it for alternate numpy backends.
                backing = getattr(value, "_mmap", None)
                if backing is not None:
                    backing.close()
                if name in locals():
                    del locals()[name]
        if parquet_tmp.exists():
            parquet_tmp.unlink()
        if not keep_work:
            for scratch in (sums_path, counts_path):
                scratch.unlink(missing_ok=True)
        if created_temp and not keep_work:
            shutil.rmtree(temporary_root, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="wide gene-by-sample TCGA RNA TSV")
    parser.add_argument("--output", required=True, help="destination patient-by-gene parquet")
    parser.add_argument("--chunk-genes", type=int, default=128, help="source gene rows parsed at once")
    parser.add_argument("--patient-chunk", type=int, default=512, help="patient rows written per parquet row group")
    parser.add_argument("--work-dir", default="", help="scratch directory for disk-backed reduction arrays")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()
    manifest = prepare_pancan_rna(
        args.source, args.output, chunk_genes=args.chunk_genes, patient_chunk=args.patient_chunk,
        work_dir=args.work_dir or None, overwrite=args.overwrite, keep_work=args.keep_work,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
