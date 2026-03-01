#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


@dataclass(frozen=True)
class LTRHit:
    index: int
    seqid: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    strand: str # '+' or '-'
    score: Optional[float] = None
    similarity: Optional[float] = None


INDEX_RE = re.compile(r"^\[\s*(\d+)\s*\]$")


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def parse_index(field0: str) -> Optional[int]:
    m = INDEX_RE.match(field0.strip())
    if not m:
        return None
    return int(m.group(1))


def parse_location(loc: str) -> Optional[Tuple[int, int]]:
    # expects "start-end"
    if "-" not in loc:
        return None
    a, b = loc.split("-", 1)
    try:
        start = int(a)
        end = int(b)
        return start, end
    except ValueError:
        return None


def normalize_strand(s: str) -> Optional[str]:
    s = s.strip()
    if s in {"+", "plus", "Plus"}:
        return "+"
    if s in {"-", "minus", "Minus"}:
        return "-"
    # Some tools use "D"/"C" or "F"/"R" in certain outputs; extend here if needed.
    return None


def parse_ltrfinder_table(path: Path) -> List[LTRHit]:
    """
    Parse LTR-Finder tabular output.

    IMPORTANT:
    LTR-Finder outputs can vary. This parser assumes:
      - Hit lines start with [n]
      - Column 2 is SeqID
      - Column 3 is 'start-end'
      - Strand/score/similarity are present somewhere later.
    If your output differs, adjust STRAND_COL / SCORE_COL / SIMILARITY_COL below.
    """
    hits: List[LTRHit] = []

    # These are the risky assumptions in your original script:
    STRAND_COL = 12
    SCORE_COL = 13
    SIMILARITY_COL = -1

    with path.open("r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line.startswith("["):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                warn(f"{path.name}:{ln}: Too few columns, skipping.")
                continue

            idx = parse_index(parts[0])
            if idx is None:
                warn(f"{path.name}:{ln}: Could not parse hit index from '{parts[0]}', skipping.")
                continue

            seqid = parts[1].strip()
            loc = parse_location(parts[2].strip())
            if loc is None:
                warn(f"{path.name}:{ln}: Could not parse location from '{parts[2]}', skipping.")
                continue
            start, end = loc

            if start <= 0 or end <= 0:
                warn(f"{path.name}:{ln}: Non-positive coordinates {start}-{end}, skipping.")
                continue

            # Pull optional fields
            strand = None
            score = None
            similarity = None

            try:
                strand = normalize_strand(parts[STRAND_COL])
            except Exception:
                strand = None

            try:
                score = float(parts[SCORE_COL])
            except Exception:
                score = None

            try:
                similarity = float(parts[SIMILARITY_COL])
            except Exception:
                similarity = None

            if strand is None:
                # If strand missing, default to '+' but warn
                warn(f"{path.name}:{ln}: Strand not found/recognized; defaulting to '+'.")
                strand = "+"

            # Ensure start <= end
            if start > end:
                start, end = end, start

            hits.append(LTRHit(idx, seqid, start, end, strand, score, similarity))

    return hits


def clamp_coords(start: int, end: int, seqlen: int) -> Tuple[int, int]:
    # start/end are 1-based inclusive
    start = max(1, min(seqlen, start))
    end = max(1, min(seqlen, end))
    if start > end:
        start, end = end, start
    return start, end


def extract_hits_from_fasta(hits: List[LTRHit], genome_fasta: Path, output_fasta: Path) -> Tuple[int, int]:
    """
    Extract hits from genome FASTA.
    """
    genome = SeqIO.index(str(genome_fasta), "fasta")
    written = 0
    skipped = 0

    with output_fasta.open("w", encoding="utf-8") as out:
        for h in hits:
            if h.seqid not in genome:
                warn(f"SeqID '{h.seqid}' not found in genome FASTA; skipping hit {h.index}.")
                skipped += 1
                continue

            rec = genome[h.seqid]
            seqlen = len(rec.seq)

            start, end = clamp_coords(h.start, h.end, seqlen)

            # Convert 1-based inclusive to Python slicing (0-based, end-exclusive)
            subseq: Seq = rec.seq[start - 1 : end]

            if h.strand == "-":
                subseq = subseq.reverse_complement()

            # Header includes optional fields if present
            extra = []
            if h.score is not None:
                extra.append(f"score={h.score}")
            if h.similarity is not None:
                extra.append(f"sim={h.similarity}")
            extra_str = ("|" + "|".join(extra)) if extra else ""

            out.write(f">{h.seqid}_{h.index}_{start}_{end}_strand{h.strand}{extra_str}\n{subseq}\n")
            written += 1

    genome.close()
    return written, skipped


def main() -> int:
    ap = argparse.ArgumentParser(description="Extract LTR-Finder hits from genome FASTA (tabular output).")
    ap.add_argument("-l", "--ltr", required=True, type=Path, help="LTR-Finder tabular output file")
    ap.add_argument("-g", "--genome", required=True, type=Path, help="Genome assembly FASTA file")
    ap.add_argument("-o", "--output", required=True, type=Path, help="Output FASTA file for extracted hits")
    args = ap.parse_args()

    if not args.ltr.exists():
        ap.error(f"LTR-Finder output not found: {args.ltr}")
    if not args.genome.exists():
        ap.error(f"Genome FASTA not found: {args.genome}")

    hits = parse_ltrfinder_table(args.ltr)
    print(f"Parsed {len(hits)} LTR-Finder hits")

    written, skipped = extract_hits_from_fasta(hits, args.genome, args.output)
    print(f"Extracted sequences written: {written}")
    if skipped:
        print(f"Skipped hits: {skipped}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())