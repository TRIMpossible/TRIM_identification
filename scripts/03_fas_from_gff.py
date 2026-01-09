#!/usr/bin/env python3
"""Extract FASTA sequences from GFF features with flanking sequence and optional splitting.

Requires:
  - biopython
  - bcbio-gff (provides BCBio.GFF)

Input:
  - genome FASTA (reference)
  - GFF files (typically merged output from step 02)
    with attributes: ID=..., middle_pos=..., optional cut_pos1/cut_pos2

Output:
  - For each flank size, two FASTA files per input GFF:
    * *_non_split.fasta  (full region including flanks)
    * *_split.fasta      (split parts: 2 parts for small regions, 4 parts for large regions if cut_pos1/2 present)

Notes on coordinates:
  - BCBio.GFF parses GFF into Biopython SeqFeatures with 0-based, end-exclusive locations.
  - Your GFF coordinates are preserved in output headers as provided by the parser (start/end integers).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from BCBio import GFF

def extract_sequence_with_flanking(sequence, start: int, end: int, strand: int, flank: int):
    start_flank = max(0, start - flank)
    end_flank = min(len(sequence), end + flank)
    seq = sequence[start_flank:end_flank]
    return seq if strand == 1 else seq.reverse_complement(), start_flank, end_flank

def split_four(sequence, start_flank: int, middle_pos: int, cut_pos1: int, cut_pos2: int, strand: int):
    # Convert genome positions to indices within `sequence`
    m = middle_pos - start_flank
    c1 = cut_pos1 - start_flank
    c2 = cut_pos2 - start_flank
    # Clamp defensively
    m = max(0, min(len(sequence), m))
    c1 = max(0, min(len(sequence), c1))
    c2 = max(0, min(len(sequence), c2))
    # Ensure order
    c1, m, c2 = sorted([c1, m, c2])

    parts = [sequence[:c1], sequence[c1:m], sequence[m:c2], sequence[c2:]]
    if strand != 1:
        parts = [p.reverse_complement() for p in parts]
    return parts

def split_two(sequence, start_flank: int, middle_pos: int, strand: int):
    m = middle_pos - start_flank
    m = max(0, min(len(sequence), m))
    parts = [sequence[:m], sequence[m:]]
    if strand != 1:
        parts = [p.reverse_complement() for p in parts]
    return parts

def main() -> int:
    ap = argparse.ArgumentParser(description="Extract FASTA from GFF features with flanks and optional splitting.")
    ap.add_argument("--input-dir", "-i", required=True, type=Path, help="Directory with input .gff files.")
    ap.add_argument("--genome", "-g", required=True, type=Path, help="Genome/reference FASTA file.")
    ap.add_argument("--outdir", "-o", required=True, type=Path, help="Output directory.")
    ap.add_argument("--pattern", default="*_merged.gff", help="Glob pattern for GFF inputs. Default: *_merged.gff")
    ap.add_argument(
        "--flanks",
        default="2000",
        help="Comma-separated flank sizes in bp (e.g. '500,2000'). Default: 2000",
    )
    ap.add_argument(
        "--split-threshold",
        type=int,
        default=10000,
        help="If feature length >= this and cut_pos1/2 exist, split into 4 parts. Default: 10000",
    )

    args = ap.parse_args()

    if not args.input_dir.exists():
        ap.error(f"Input dir not found: {args.input_dir}")
    if not args.genome.exists():
        ap.error(f"Genome FASTA not found: {args.genome}")

    flank_sizes: List[int] = []
    for part in str(args.flanks).split(","):
        part = part.strip()
        if not part:
            continue
        flank_sizes.append(int(part))
    if not flank_sizes:
        ap.error("No valid flank sizes provided.")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load genome into dict once
    genome_dict = SeqIO.to_dict(SeqIO.parse(str(args.genome), "fasta"))

    gff_files = sorted(args.input_dir.glob(args.pattern))
    if not gff_files:
        ap.error(f"No GFF files matched '{args.pattern}' in {args.input_dir}")

    for gff_file in gff_files:
        base = gff_file.stem
        print(f"Processing: {gff_file.name}")

        with gff_file.open("r", encoding="utf-8") as gff_handle:
            for rec in GFF.parse(gff_handle, base_dict=genome_dict):
                for feature in rec.features:
                    feature_id = feature.qualifiers.get("ID", ["unknown"])[0]
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = int(feature.location.strand or 1)

                    middle_pos = feature.qualifiers.get("middle_pos", [None])[0]
                    if middle_pos is None:
                        # Without a middle_pos, splitting isn't well-defined; still write non-split.
                        middle_pos_int = start + (end - start) // 2
                    else:
                        middle_pos_int = int(middle_pos)

                    cut_pos1 = feature.qualifiers.get("cut_pos1", [None])[0]
                    cut_pos2 = feature.qualifiers.get("cut_pos2", [None])[0]
                    cut_pos1_int = int(cut_pos1) if cut_pos1 is not None else None
                    cut_pos2_int = int(cut_pos2) if cut_pos2 is not None else None

                    length = end - start

                    for flank in flank_sizes:
                        flank_dir = args.outdir / f"flank{flank}"
                        flank_dir.mkdir(parents=True, exist_ok=True)

                        non_split_path = flank_dir / f"{base}_flank{flank}_non_split.fasta"
                        split_path = flank_dir / f"{base}_flank{flank}_split.fasta"

                        full_seq, start_flank, end_flank = extract_sequence_with_flanking(rec.seq, start, end, strand, flank)

                        with non_split_path.open("a", encoding="utf-8") as h:
                            h.write(f">{feature_id}|{rec.id}:{start}-{end}|strand{strand}|flank{flank}\n{full_seq}\n")

                        # Splitting
                        if length >= args.split_threshold and cut_pos1_int is not None and cut_pos2_int is not None:
                            parts = split_four(full_seq, start_flank, middle_pos_int, cut_pos1_int, cut_pos2_int, strand)
                        else:
                            parts = split_two(full_seq, start_flank, middle_pos_int, strand)

                        with split_path.open("a", encoding="utf-8") as h:
                            for idx, part_seq in enumerate(parts, start=1):
                                if len(part_seq) == 0:
                                    continue
                                h.write(f">{feature_id}|{rec.id}:{start}-{end}|strand{strand}|flank{flank}|part{idx}\n{part_seq}\n")

    print("Done.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
