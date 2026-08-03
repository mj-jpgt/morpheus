"""Choose which TCGA solid-tissue-normal WSIs to patch for the dilution curve.

Every ``.svs`` under ``tcga_normal_recovery/`` and ``tcga_normal_modalities/`` in
``mj0jpgg/tcga_cptac_cgga`` carries TCGA sample-type code ``11`` (solid tissue
normal).  The study of origin is not in the filename, so it is recovered from the
tissue-source-site field of the barcode via a majority map built from the frozen
cancer registry.  Slides whose TSS is not in the registry are kept but marked
``cancer=UNKNOWN`` -- they still belong in the pan-cancer pool.

Two pools are needed:
  * a pan-cancer pool, used by every patient in the primary arm;
  * per-cancer pools for the organ-matched sensitivity arm, which only exist for
    the test cancers that have enough normal donors.
"""
from __future__ import annotations

import argparse

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", required=True, help="tcga_cancer_registry.parquet")
    parser.add_argument("--output", required=True)
    parser.add_argument("--matched-cancers", default="")
    parser.add_argument("--per-matched-cancer", type=int, default=30)
    parser.add_argument("--per-other-cancer", type=int, default=8)
    parser.add_argument("--max-size-mb", type=float, default=500.0)
    args = parser.parse_args()

    from huggingface_hub import HfApi, get_token

    info = HfApi(token=get_token()).repo_info("mj0jpgg/tcga_cptac_cgga", repo_type="dataset",
                                              files_metadata=True)
    slides = [(s.rfilename, s.size) for s in info.siblings
              if s.rfilename.endswith(".svs") and s.rfilename.startswith("tcga_normal_")]

    registry = pd.read_parquet(args.registry)
    registry["tss"] = registry.patient_id.str.split("-").str[1]
    tss_map = {tss: group.cancer_type.value_counts().index[0] for tss, group in registry.groupby("tss")}

    rows = []
    for name, size in slides:
        barcode = name.split("/")[-1].split(".")[0]
        pieces = barcode.split("-")
        if len(pieces) < 4 or pieces[3][:2] != "11":
            continue
        rows.append({"file": name, "size": int(size or 0), "patient": "-".join(pieces[:3]),
                     "tss": pieces[1], "cancer": tss_map.get(pieces[1], "UNKNOWN"), "slide": barcode})
    frame = pd.DataFrame(rows).drop_duplicates("slide")
    frame = frame[frame["size"] <= args.max_size_mb * 1e6]

    matched = {c for c in args.matched_cancers.split(",") if c}
    keep = []
    for cancer, group in frame.groupby("cancer"):
        budget = args.per_matched_cancer if cancer in matched else args.per_other_cancer
        # Smallest-first: a normal slide's tissue content does not scale with its
        # file size, so this buys donor breadth per byte transferred.
        keep.append(group.sort_values("size").head(budget))
    selected = pd.concat(keep).sort_values(["cancer", "size"]).reset_index(drop=True)
    selected.to_csv(args.output, index=False)
    print(selected.cancer.value_counts().to_string())
    print(f"[selected] {len(selected)} slides, {selected['size'].sum() / 1e9:.1f} GB, "
          f"{selected.patient.nunique()} donors -> {args.output}")


if __name__ == "__main__":
    main()
