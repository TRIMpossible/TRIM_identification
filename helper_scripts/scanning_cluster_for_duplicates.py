#!/usr/bin/env python3
"""
scanning_cluster_for_duplicates.py

Aligns each FASTA file in an input directory using MAFFT, computes pairwise
identity, detects overlapping genomic coordinates (based on sequence IDs),
writes logs, and optionally writes a reduced alignment with overlaps removed.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple, Set, Dict

from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord


# ----------------------------
# CLI / config
# ----------------------------

@dataclass
class Config:
    input_dir: Path
    output_dir: Path
    mafft: str
    mafft_args: List[str]
    min_identity: Optional[float]
    overlap_policy: str
    extensions: Tuple[str, ...]


def parse_args() -> Config:
    ap = argparse.ArgumentParser(
        description="Align cluster FASTAs with MAFFT, compute pairwise identity, detect coordinate overlaps, write reduced alignments."
    )
    ap.add_argument("--input-dir", "-i", required=True, type=Path, help="Directory containing cluster FASTA files.")
    ap.add_argument("--output-dir", "-o", required=True, type=Path, help="Directory for alignments/logs/reduced outputs.")
    ap.add_argument("--mafft", default=None, help="Path to MAFFT executable. If omitted, uses MAFFT in PATH.")
    ap.add_argument(
        "--mafft-args",
        default="--auto",
        help='Extra MAFFT args as a single string (default: "--auto"). Example: "--auto --thread 4"',
    )
    ap.add_argument(
        "--min-identity",
        type=float,
        default=None,
        help="Optional: flag pairs below this percent identity in the log (does not remove sequences).",
    )
    ap.add_argument(
        "--overlap-policy",
        choices=["keep-first", "keep-longest"],
        default="keep-first",
        help="How to resolve overlaps in reduced output: keep-first (default) or keep-longest.",
    )
    ap.add_argument(
        "--ext",
        default=".fasta,.fa,.fna",
        help='Comma-separated FASTA extensions to process (default: ".fasta,.fa,.fna").',
    )

    args = ap.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    mafft = args.mafft or shutil.which("mafft")

    if not input_dir.exists():
        ap.error(f"Input directory not found: {input_dir}")
    if mafft is None:
        ap.error("MAFFT not found. Install it or provide --mafft /path/to/mafft")

    exts = tuple(e.strip().lower() for e in args.ext.split(",") if e.strip())
    mafft_args = args.mafft_args.strip().split() if args.mafft_args.strip() else []

    return Config(
        input_dir=input_dir,
        output_dir=output_dir,
        mafft=str(mafft),
        mafft_args=mafft_args,
        min_identity=args.min_identity,
        overlap_policy=args.overlap_policy,
        extensions=exts,
    )


# ----------------------------
# Coordinate parsing + overlap
# ----------------------------

_COORD_RE = re.compile(
    # tries to match ..._<chrom>_<start>_<end>... anywhere in the ID
    # e.g. Fam1_00519_CANROG010000003.1_98624_99035_strand-1_flank0_1
    r"(?P<chrom>[^_]+)_(?P<start>\d+)_(?P<end>\d+)"
)

def extract_coords(seq_id: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """
    Extract (chrom, start, end) from a sequence ID.

    Supports the common pattern: <chrom>_<start>_<end> somewhere in the ID.
    Returns (None, None, None) if not found.
    """
    m = _COORD_RE.search(seq_id)
    if not m:
        return None, None, None
    chrom = m.group("chrom")
    start = int(m.group("start"))
    end = int(m.group("end"))
    # normalize order
    if start > end:
        start, end = end, start
    return chrom, start, end

def overlaps(chrom1: str, s1: int, e1: int, chrom2: str, s2: int, e2: int) -> bool:
    return chrom1 == chrom2 and s1 <= e2 and s2 <= e1

def span_len(s: int, e: int) -> int:
    return abs(e - s) + 1


# ----------------------------
# Identity calculation
# ----------------------------

def calculate_pairwise_identity(alignment_path: Path) -> Tuple[List[Tuple[str, str, float]], float]:
    """
    Computes pairwise % identity ignoring gap positions (only positions where both residues are not '-').
    Returns (pairwise list, average identity).
    """
    alignment = AlignIO.read(str(alignment_path), "fasta")
    results: List[Tuple[str, str, float]] = []
    total = 0.0
    n = 0

    for a, b in combinations(alignment, 2):
        matches = 0
        length = 0
        for ra, rb in zip(str(a.seq), str(b.seq)):
            if ra == "-" or rb == "-":
                continue
            length += 1
            if ra == rb:
                matches += 1
        ident = (matches / length) * 100.0 if length else 0.0
        results.append((a.id, b.id, ident))
        total += ident
        n += 1

    avg = total / n if n else 0.0
    return results, avg


# ----------------------------
# MAFFT runner
# ----------------------------

def run_mafft(mafft: str, mafft_args: Sequence[str], input_fasta: Path, output_fasta: Path) -> None:
    """
    Run MAFFT and write aligned FASTA.
    """
    cmd = [mafft, *mafft_args, str(input_fasta)]
    # Write stdout directly to file
    with output_fasta.open("w", encoding="utf-8") as out:
        cp = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if cp.returncode != 0:
        raise RuntimeError(f"MAFFT failed ({input_fasta.name}):\n{cp.stderr}")


# ----------------------------
# Overlap resolution
# ----------------------------

def pick_overlap_keeps(
    records: List[SeqRecord],
    policy: str,
) -> Set[str]:
    """
    Returns a set of IDs to drop due to overlaps.
    Overlaps are detected based on coords parsed from record.id.
    """
    coords: List[Tuple[SeqRecord, str, int, int]] = []
    for r in records:
        chrom, s, e = extract_coords(r.id)
        if chrom is None:
            continue
        coords.append((r, chrom, s, e))

    drop: Set[str] = set()
    # Compare all pairs with coordinates
    for i in range(len(coords)):
        r1, c1, s1, e1 = coords[i]
        if r1.id in drop:
            continue
        for j in range(i + 1, len(coords)):
            r2, c2, s2, e2 = coords[j]
            if r2.id in drop:
                continue
            if overlaps(c1, s1, e1, c2, s2, e2):
                # Decide which to drop
                if policy == "keep-first":
                    drop.add(r2.id)
                elif policy == "keep-longest":
                    len1 = span_len(s1, e1)
                    len2 = span_len(s2, e2)
                    # keep longer span; tie -> keep first
                    if len2 > len1:
                        drop.add(r1.id)
                        break
                    else:
                        drop.add(r2.id)
    return drop


# ----------------------------
# Main processing
# ----------------------------

def main() -> int:
    cfg = parse_args()

    align_dir = cfg.output_dir
    log_dir = cfg.output_dir / "logs"
    reduced_dir = cfg.output_dir / "reduced"

    align_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    reduced_dir.mkdir(parents=True, exist_ok=True)

    print(f"MAFFT: {cfg.mafft}")
    print(f"Input: {cfg.input_dir}")
    print(f"Output: {cfg.output_dir}")
    print(f"Logs: {log_dir}")
    print(f"Reduced: {reduced_dir}")

    fasta_files = sorted(
        p for p in cfg.input_dir.iterdir()
        if p.is_file() and p.suffix.lower() in cfg.extensions
    )

    if not fasta_files:
        print("No FASTA files found in input directory.")
        return 0

    for fasta in fasta_files:
        base = fasta.stem
        aligned_path = align_dir / f"{base}.aligned.fasta"
        log_path = log_dir / f"{base}.log"
        tsv_path = log_dir / f"{base}.pairwise_identity.tsv"
        reduced_path = reduced_dir / f"{base}.reduced.fasta"

        print(f"\nProcessing: {fasta.name}")

        try:
            run_mafft(cfg.mafft, cfg.mafft_args, fasta, aligned_path)
            print(f"✅ Alignment written: {aligned_path}")

            records = list(SeqIO.parse(str(aligned_path), "fasta"))
            pair_results, avg_identity = calculate_pairwise_identity(aligned_path)

            # Overlap detection + reduction
            drop_ids = pick_overlap_keeps(records, policy=cfg.overlap_policy)
            kept = [r for r in records if r.id not in drop_ids]
            SeqIO.write(kept, str(reduced_path), "fasta")

            # Write logs
            with log_path.open("w", encoding="utf-8") as log:
                log.write(f"Input: {fasta}\n")
                log.write(f"Aligned: {aligned_path}\n")
                log.write(f"Policy: {cfg.overlap_policy}\n")
                log.write(f"Sequences in alignment: {len(records)}\n")
                log.write(f"Sequences kept (reduced): {len(kept)}\n")
                log.write(f"Sequences removed: {len(drop_ids)}\n\n")
                log.write(f"Average pairwise identity: {avg_identity:.2f}%\n\n")

                # Identity warnings
                if cfg.min_identity is not None:
                    low = [(a, b, pid) for a, b, pid in pair_results if pid < cfg.min_identity]
                    log.write(f"Pairs below {cfg.min_identity:.2f}% identity: {len(low)}\n")
                    for a, b, pid in low[:200]:
                        log.write(f"  LOW_ID\t{a}\t{b}\t{pid:.2f}\n")
                    if len(low) > 200:
                        log.write("  ... (truncated)\n")
                    log.write("\n")

                # Overlap details
                if drop_ids:
                    log.write("Removed due to coordinate overlaps:\n")
                    for did in sorted(drop_ids):
                        chrom, s, e = extract_coords(did)
                        log.write(f"  {did}\t{chrom}:{s}-{e}\n")
                    log.write("\n")

            with tsv_path.open("w", encoding="utf-8") as tsv:
                tsv.write("seq1\tseq2\tidentity\n")
                for a, b, pid in pair_results:
                    tsv.write(f"{a}\t{b}\t{pid:.4f}\n")

            print(f"📝 Log written: {log_path}")
            print(f"🧾 Pairwise identity TSV: {tsv_path}")
            print(f"📦 Reduced alignment written: {reduced_path}")

        except Exception as e:
            print(f"❌ Error processing {fasta.name}: {e}")

    print("\n✅ Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())