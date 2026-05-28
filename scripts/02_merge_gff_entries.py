#!/usr/bin/env python3
"""Merge nearby GFF features and keep attributes from the largest hit.

Reads GFF3 files, filters by identity threshold (identity= in attributes),
merges features that are within --merge-distance, and writes *_merged.gff
files plus per-input log files.

Assumptions
-----------
- Input features have attributes including:
  - ID=...
  - identity=<float>  (percentage)
- Coordinates are treated as 1-based inclusive (standard GFF). We only
  compare distances between features, so indexing convention does not
  affect merging behavior.
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, List, Optional, Tuple

Feature = Tuple[int, int, str, str, str, str]
# start, end, feature_type, strand, attributes, feature_id

def extract_identity(attributes: str) -> Optional[float]:
    for attribute in attributes.split(";"):
        if attribute.startswith("identity="):
            try:
                return float(attribute.split("=", 1)[1])
            except ValueError:
                return None
    return None

def extract_id(attributes: str) -> Optional[str]:
    for attribute in attributes.split(";"):
        if attribute.startswith("ID="):
            return attribute.split("=", 1)[1]
    return None

def update_id_and_middle(attributes: str, new_id: str, middle_pos: Optional[int]) -> str:
    parts = attributes.split(";") if attributes else []
    out = []
    seen_id = False
    seen_middle = False
  
    for part in parts:
        if part.startswith("ID="):
            out.append(f"ID={new_id}")
            seen_id = True
        elif part.startswith("middle_pos="):
            if middle_pos is not None:
                out.append(f"middle_pos={middle_pos}")
            seen_middle = True
        else:
            out.append(part)

    if not seen_id:
        out.insert(0, f"ID={new_id}")
    if middle_pos is not None and not seen_middle:
        out.append(f"middle_pos={middle_pos}")

    # Remove any empty fields
    out = [p for p in out if p]
    return ";".join(out)

def calculate_cut_positions(start: int, end: int, split_threshold: int) -> List[int]:
    length = end - start
    if length >= split_threshold:
        return [start + length // 4 * i for i in range(1, 4)]
    return []

def merge_gff_entries(
    gff_file: Path,
    output_file: Path,
    merge_distance: int,
    identity_threshold: float,
    split_threshold: int,
    log_dir: Path,
) -> None:
    entries: DefaultDict[Tuple[str, str], List[Feature]] = defaultdict(list)
    merged_features_dict: DefaultDict[Tuple[str, str], List[Tuple[int, int, str, str, str, List[str]]]] = defaultdict(list)

    # Read features
    with gff_file.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts[:9]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue

            if end_i < start_i:
                start_i, end_i = end_i, start_i
          
            identity = extract_identity(attributes)
            feature_id = extract_id(attributes) or "unknown"

            if identity is None or identity < identity_threshold:
                continue

            entries[seq_id, strand].append((start_i, end_i, feature_type, strand, attributes, feature_id))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8") as out:
        out.write("##gff-version 3\n")
        
        for (seq_id, strand_key), features in sorted(entries.items()):
            if not features:
                continue
            features.sort(key=lambda x: x[0])

            current_start, current_end, current_type, current_strand, largest_hit_attributes, current_id = features[0]
            largest_hit_size = current_end - current_start
            merged_ids: List[str] = [current_id]

            def flush_group() -> None:
                merged_middle_pos = current_start + (current_end - current_start) // 2
                updated_attributes = update_id_and_middle(largest_hit_attributes, current_id, merged_middle_pos)

                cut_positions = calculate_cut_positions(current_start, current_end, split_threshold)
                cut_positions = [pos for pos in cut_positions if pos != merged_middle_pos]
              
                if cut_positions:
                    cut_pos_str = ";".join(f"cut_pos{i}={pos}" for i, pos in enumerate(cut_positions, start=1))
                    updated_attributes = updated_attributes + ";" + cut_pos_str

                out.write(f"{seq_id}\t.\t{current_type}\t{current_start}\t{current_end}\t.\t"
                          f"{current_strand}\t.\t{updated_attributes}\n"
                merged_features_dict[(seq_id, current_strand)].append((current_start, current_end, current_type, current_strand, updated_attributes, merged_ids.copy()))

            for start, end, feature_type, strand, attributes, feature_id in features[1:]:
                current_hit_size = end - start

                same_strand = strand ==current_strand
                close_enough = start -current_end <= merge_distance
              
                if current_hit_size > largest_hit_size:
                    largest_hit_size = current_hit_size
                    largest_hit_attributes = attributes
                    current_id = feature_id

             else:
                flush_group()
                current_start, current_end = start, end
                current_type, current_strand = feature_type, strand
                largest_hit_attributes = attributes
                current_id = feature_id
                largest_hit_size = current_hit_size
                merged_ids = [feature_id]

            flush_group()

    write_log(entries, merged_features_dict, gff_file.name, log_dir, identity_threshold, merge_distance)
    print(f"Merged: {gff_file.name} -> {output_file.name}")

def write_log(
    entries: DefaultDict[Tuple [str, str], List[Feature]],
    merged_features_dict: DefaultDict[Tuple[str, str], List[Tuple[int,int,str,str,str,List[str]]]],
    filename: str,
    log_dir: Path,
    identity_threshold: float,
    merge_distance: int,
) -> None:
    log_file = log_dir / f"{Path(filename).stem}_log.txt"

    with log_file.open("w", encoding="utf-8") as log:
        log.write(f"Input file: {filename}\n")
        log.write(f"identity_threshold: {identity_threshold}\n")
        log.write(f"merge_distance: {merge_distance}\n\n")
        log.write("Merging rule: same esequence ID, same strand, within merge distance\n\n")

        for key in sorted(entries):
            seq_id, strand = key
            original_count = len(entrie[key])
            merged_features = merged_features_dict.get(key, [])
            merged_count = len(merged_features)
            hits_merged = original_count - merged_count

            log.write(f"{seq_id} strand={strand}\n")
            for mf in merged_features:
                merged_ids = mf[5]
                merged_ids_str = ", ".join([m for mf in merged_ids if m])
                log.write(f"- {merged_ids_str}\n")
            log.write("\n")

def main() -> int:
    ap = argparse.ArgumentParser(description="Merge nearby GFF entries with identity filtering.")
    ap.add_argument("--input-dir", "-i", required=True, type=Path, help="Folder containing .gff files.")
    ap.add_argument("--output-dir", "-o", required=True, type=Path, help="Folder for *_merged.gff outputs.")
    ap.add_argument("--log-dir", required=True, type=Path, help="Folder for log files.")
    ap.add_argument("--merge-distance", type=int, default=5000, help="Max gap to merge (bp). Default: 5000.")
    ap.add_argument("--identity-threshold", type=float, default=75.0, help="Minimum identity to keep. Default: 75.")
    ap.add_argument("--split-threshold", type=int, default=10000, help="Add cut_pos* if merged span >= this. Default: 10000.")
    ap.add_argument("--pattern", default="*.gff", help="Glob pattern for input files. Default: *.gff")

    args = ap.parse_args()

    if not args.input_dir.exists():
        ap.error(f"Input dir not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(args.input_dir.glob(args.pattern))
    if not files:
        ap.error(f"No files matched '{args.pattern}' in {args.input_dir}")

    for gff in files:
        out = args.output_dir / f"{gff.stem}_merged.gff"
        merge_gff_entries(
            gff_file=gff,
            output_file=out,
            merge_distance=args.merge_distance,
            identity_threshold=args.identity_threshold,
            split_threshold=args.split_threshold,
            log_dir=args.log_dir,
        )

    print("Done.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
