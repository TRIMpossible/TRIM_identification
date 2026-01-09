#!/usr/bin/env python3
"""Split or filter BLAST tabular output by query ID.

Default behavior: read a BLAST outfmt 6/7-like tabular file (TSV),
group lines by the query column, and write one file per query.

This makes downstream steps easier because each query gets its own file.

Examples
--------
# Split all queries into separate files
python 00_filter_BLAST_by_query.py --input results.m6 --outdir filtered

# Only keep queries listed in a text file (one query ID per line)
python 00_filter_BLAST_by_query.py --input results.m6 --outdir filtered --queries keep.txt
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

def _safe_filename(name: str, max_len: int = 150) -> str:
    # Replace characters that are risky on Windows/macOS/Linux filesystems
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    safe = safe.strip("._-") or "query"
    return safe[:max_len]

def read_queries_file(path: Path) -> set[str]:
    queries: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        queries.add(line)
    return queries

def iter_blast_lines(input_file: Path) -> Iterable[str]:
    with input_file.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            # Ignore comment lines (outfmt 7 style)
            if line.startswith("#"):
                continue
            yield line

def split_blast_by_query(
    input_file: Path,
    outdir: Path,
    query_col: int = 0,
    keep_queries: Optional[set[str]] = None,
    suffix: str = "_filtered_output.txt",
) -> Dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)

    handles: Dict[str, object] = {}
    outpaths: Dict[str, Path] = {}

    try:
        for line in iter_blast_lines(input_file):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= query_col:
                continue
            query_id = fields[query_col]
            if keep_queries is not None and query_id not in keep_queries:
                continue

            if query_id not in handles:
                fname = _safe_filename(query_id) + suffix
                outpath = outdir / fname
                handles[query_id] = outpath.open("w", encoding="utf-8")
                outpaths[query_id] = outpath

            handles[query_id].write(line)
    finally:
        for h in handles.values():
            try:
                h.close()
            except Exception:
                pass

    return outpaths

def main() -> int:
    ap = argparse.ArgumentParser(description="Split/filter BLAST tabular output by query ID.")
    ap.add_argument("--input", "-i", required=True, type=Path, help="BLAST tabular output file (outfmt 6/7).")
    ap.add_argument("--outdir", "-o", required=True, type=Path, help="Output directory for per-query files.")
    ap.add_argument("--query-col", type=int, default=0, help="0-based column index for query ID. Default: 0.")
    ap.add_argument(
        "--queries",
        type=Path,
        default=None,
        help="Optional file with query IDs to keep (one per line). If omitted, all queries are written.",
    )
    ap.add_argument(
        "--suffix",
        default="_filtered_output.txt",
        help="Suffix for output files. Default: _filtered_output.txt",
    )

    args = ap.parse_args()

    if not args.input.exists():
        ap.error(f"Input file not found: {args.input}")

    keep = read_queries_file(args.queries) if args.queries else None
    outpaths = split_blast_by_query(
        input_file=args.input,
        outdir=args.outdir,
        query_col=args.query_col,
        keep_queries=keep,
        suffix=args.suffix,
    )

    print(f"Wrote {len(outpaths)} file(s) to: {args.outdir}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
