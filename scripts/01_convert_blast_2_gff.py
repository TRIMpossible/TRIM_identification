#!/usr/bin/env python3
"""Convert BLAST tabular output files (outfmt 6) into GFF3.

This script scans an input folder for BLAST tabular files and writes one
GFF per input file.

Notes
-----
- Assumes BLAST outfmt 6 with at least 12 columns:
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
- Uses sseqid as the GFF seqid.
- Determines strand by comparing sstart and send.
"""

from __future__ import annotations

import argparse
from pathlib import Path

def blast_file_to_gff(
    blast_output_file: Path,
    gff_output_file: Path,
    family_tag: str,
    feature_type: str = "alignment",
    source: str = "BLAST",
) -> int:
    serial_number = 1
    written = 0

    gff_output_file.parent.mkdir(parents=True, exist_ok=True)
    with blast_output_file.open("r", encoding="utf-8") as blast_file, gff_output_file.open("w", encoding="utf-8") as gff_file:
        gff_file.write("##gff-version 3\n")
        for line in blast_file:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            serial_id = f"{family_tag}_{serial_number:05d}"
            serial_number += 1

            subject_id = fields[1]
            try:
                identity = float(fields[2])
                start = int(fields[8])
                end = int(fields[9])
                bitscore = float(fields[11])
            except ValueError:
                continue

            if start < end:
                strand = "+"
            else:
                strand = "-"
                start, end = end, start

            attributes = f"ID={serial_id};bitscore={bitscore};identity={identity}"
            gff_entry = "\t".join([subject_id, source, feature_type, str(start), str(end), str(bitscore), strand, ".", attributes])
            gff_file.write(gff_entry + "\n")
            written += 1

    return written

def main() -> int:
    ap = argparse.ArgumentParser(description="Convert BLAST outfmt 6 files in a folder to GFF3.")
    ap.add_argument("--input-dir", "-i", required=True, type=Path, help="Folder with BLAST tabular files.")
    ap.add_argument("--output-dir", "-o", required=True, type=Path, help="Folder for generated GFF files.")
    ap.add_argument(
        "--pattern",
        default="*.txt",
        help="Glob pattern for BLAST files (e.g. '*.txt' or '*.m6'). Default: *.txt",
    )
    ap.add_argument(
        "--family-tag-from",
        choices=["prefix", "filename"],
        default="prefix",
        help="How to derive family_tag: 'prefix' uses text before first '_' in filename; 'filename' uses full stem.",
    )

    args = ap.parse_args()

    if not args.input_dir.exists():
        ap.error(f"Input dir not found: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(args.input_dir.glob(args.pattern))
    if not files:
        ap.error(f"No files matched pattern '{args.pattern}' in {args.input_dir}")

    total_written = 0
    for f in files:
        stem = f.stem
        if args.family_tag_from == "prefix":
            family_tag = stem.split("_")[0] if "_" in stem else stem
        else:
            family_tag = stem

        out = args.output_dir / (stem + ".gff")
        written = blast_file_to_gff(f, out, family_tag=family_tag)
        total_written += written
        print(f"{f.name} -> {out.name} ({written} features)")

    print(f"Done. Total features written: {total_written}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
