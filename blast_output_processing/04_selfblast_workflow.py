#!/usr/bin/env python3
"""
Dependencies
------------
- Python 3.8+
- biopython (pip install biopython)
- NCBI BLAST+ executables in PATH: blastn, makeblastdb
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# ----------------------------
# Utilities
# ----------------------------

def which_or_fail(prog: str) -> str:
    p = shutil.which(prog)
    if not p:
        raise SystemExit(f"ERROR: Required program '{prog}' not found in PATH.")
    return p


def run_cmd(cmd: Sequence[str], *, capture_stdout: bool = False) -> List[str]:
    """Run a command safely. If capture_stdout=True, return stdout lines."""
    # Show a readable command line
    print("\n>>", " ".join([shlex_quote(c) for c in cmd]))
    try:
        if capture_stdout:
            cp = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if cp.stderr.strip():
                # BLAST can be chatty; keep it but don't fail
                print(cp.stderr.strip(), file=sys.stderr)
            return [line for line in cp.stdout.splitlines() if line.strip()]
        else:
            cp = subprocess.run(cmd, check=True)
            return []
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"ERROR: Command failed with exit code {e.returncode}: {' '.join(cmd)}")


def shlex_quote(s: str) -> str:
    """Minimal shell quoting for display only (not used for execution)."""
    if re.search(r"[ \t\n\"'\\$&;()<>|]", s):
        return '"' + s.replace('"', '\\"') + '"'
    return s


def list_fastas(d: Path) -> List[Path]:
    exts = {".fa", ".fasta", ".fna"}
    files = [p for p in d.iterdir() if p.is_file() and p.suffix.lower() in exts]
    return sorted(files)


def fasta_id(rec: SeqRecord) -> str:
    # Use "record.id" (first token up to whitespace), which is what BLAST will use.
    return rec.id


def extract_group_id(seq_id: str) -> str:
    # Remove trailing _<digits> to form a group key, else keep original.
    return re.sub(r"_[0-9]+$", "", seq_id)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


# ----------------------------
# BLAST + parsing
# ----------------------------

BLAST_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
)

@dataclass(frozen=True)
class BlastHit:
    qid: str
    sid: str
    pident: str
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: str

    @staticmethod
    def from_outfmt6_line(line: str) -> "BlastHit":
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            raise ValueError("BLAST line has <12 columns")
        return BlastHit(
            qid=f[0],
            sid=f[1],
            pident=f[2],
            length=int(f[3]),
            mismatch=int(f[4]),
            gapopen=int(f[5]),
            qstart=int(f[6]),
            qend=int(f[7]),
            sstart=int(f[8]),
            send=int(f[9]),
            evalue=f[10],
            bitscore=f[11],
        )


def blast_pair(blastn: str, query_fa: Path, subject_fa: Path) -> List[str]:
    cmd = [
        blastn,
        "-query", str(query_fa),
        "-subject", str(subject_fa),
        "-outfmt", f"6 {BLAST_FIELDS}",
        "-parse_deflines",
    ]
    return run_cmd(cmd, capture_stdout=True)


# ----------------------------
# GFF writing
# ----------------------------

def gff_line(
    seqid: str,
    source: str,
    ftype: str,
    start: int,
    end: int,
    score: str,
    strand: str,
    phase: str,
    attrs: str,
) -> str:
    return "\t".join([seqid, source, ftype, str(start), str(end), str(score), strand, phase, attrs])


def blast_hit_to_gff_per_sequence(hit: BlastHit) -> List[str]:
    # Determine strand by coordinate direction (BLAST coordinates are 1-based inclusive)
    qstrand = "+" if hit.qstart <= hit.qend else "-"
    sstrand = "+" if hit.sstart <= hit.send else "-"

    qs, qe = (hit.qstart, hit.qend) if hit.qstart <= hit.qend else (hit.qend, hit.qstart)
    ss, se = (hit.sstart, hit.send) if hit.sstart <= hit.send else (hit.send, hit.sstart)

    qattrs = f"ID={hit.qid}_vs_{hit.sid};Target={hit.sid};pident={hit.pident};evalue={hit.evalue}"
    sattrs = f"ID={hit.sid}_vs_{hit.qid};Target={hit.qid};pident={hit.pident};evalue={hit.evalue}"

    return [
        gff_line(hit.qid, "BLAST", "match", qs, qe, hit.bitscore, qstrand, ".", qattrs),
        gff_line(hit.sid, "BLAST", "match", ss, se, hit.bitscore, sstrand, ".", sattrs),
    ]


def blast_hit_to_gff_combined(hit: BlastHit, query_len: int, subject_len: int) -> List[str]:
    # Combined coordinate space: query [1..query_len], subject [query_len+1 .. query_len+subject_len]
    qstrand = "+" if hit.qstart <= hit.qend else "-"
    sstrand = "+" if hit.sstart <= hit.send else "-"

    qs, qe = (hit.qstart, hit.qend) if hit.qstart <= hit.qend else (hit.qend, hit.qstart)
    ss, se = (hit.sstart, hit.send) if hit.sstart <= hit.send else (hit.send, hit.sstart)

    ss2 = query_len + ss
    se2 = query_len + se

    seqid = f"{hit.qid}__vs__{hit.sid}"

    qattrs = f"ID={hit.qid}_hit;part=query;Target={hit.sid};pident={hit.pident};evalue={hit.evalue}"
    sattrs = f"ID={hit.sid}_hit;part=subject;Target={hit.qid};pident={hit.pident};evalue={hit.evalue}"

    return [
        gff_line(seqid, "BLAST", "query_match", qs, qe, hit.bitscore, qstrand, ".", qattrs),
        gff_line(seqid, "BLAST", "subject_match", ss2, se2, hit.bitscore, sstrand, ".", sattrs),
    ]


# ----------------------------
# Core pipeline steps
# ----------------------------

def load_fasta_records(fasta_file: Path) -> List[SeqRecord]:
    return list(SeqIO.parse(str(fasta_file), "fasta"))


def write_single_record_fasta(rec: SeqRecord, path: Path) -> None:
    with path.open("w", encoding="utf-8") as h:
        SeqIO.write(rec, h, "fasta")


def process_one_fasta(
    fasta_file: Path,
    out_dir: Path,
    blastn: str,
    write_combined_gff: bool = True,
) -> None:
    base = fasta_file.stem
    out_blast = out_dir / f"{base}_selfblast.tsv"
    out_gff = out_dir / f"{base}_selfblast.gff"
    out_gff_comb = out_dir / f"{base}_selfblast_combined.gff"

    records = load_fasta_records(fasta_file)
    if len(records) < 2:
        print(f"Skipping {fasta_file.name} (less than 2 sequences).")
        return

    # Index records by ID for lengths and temp writing
    rec_by_id: Dict[str, SeqRecord] = {fasta_id(r): r for r in records}

    # Group ids
    group_map: Dict[str, List[str]] = {}
    for rid in rec_by_id.keys():
        gid = extract_group_id(rid)
        group_map.setdefault(gid, []).append(rid)

    all_blast_lines: List[str] = []
    all_gff_lines: List[str] = []
    all_gff_comb_lines: List[str] = []

    with tempfile.TemporaryDirectory(prefix=f"selfblast_{base}_") as tmpdir:
        tmpdir_p = Path(tmpdir)

        for gid, ids in group_map.items():
            if len(ids) < 2:
                continue

            # Write each record to its own FASTA in tmp
            fa_paths: Dict[str, Path] = {}
            for rid in ids:
                p = tmpdir_p / f"{rid}.fasta"
                write_single_record_fasta(rec_by_id[rid], p)
                fa_paths[rid] = p

            # Pairwise combinations
            for qid, sid in itertools.combinations(ids, 2):
                blast_lines = blast_pair(blastn, fa_paths[qid], fa_paths[sid])
                if not blast_lines:
                    continue

                all_blast_lines.extend(blast_lines)

                for line in blast_lines:
                    try:
                        hit = BlastHit.from_outfmt6_line(line)
                    except Exception:
                        continue

                    all_gff_lines.extend(blast_hit_to_gff_per_sequence(hit))

                    if write_combined_gff:
                        qlen = len(rec_by_id[qid].seq)
                        slen = len(rec_by_id[sid].seq)
                        all_gff_comb_lines.extend(blast_hit_to_gff_combined(hit, qlen, slen))

    # Write outputs
    out_blast.write_text("\n".join(all_blast_lines) + ("\n" if all_blast_lines else ""), encoding="utf-8")
    out_gff.write_text("##gff-version 3\n" + "\n".join(all_gff_lines) + ("\n" if all_gff_lines else ""), encoding="utf-8")

    if write_combined_gff:
        out_gff_comb.write_text("##gff-version 3\n" + "\n".join(all_gff_comb_lines) + ("\n" if all_gff_comb_lines else ""), encoding="utf-8")

    print(f"Wrote: {out_blast}")
    print(f"Wrote: {out_gff}")
    if write_combined_gff:
        print(f"Wrote: {out_gff_comb}")


def read_gff_seqids(gff_file: Path) -> Set[str]:
    seqids: Set[str] = set()
    with gff_file.open("r", encoding="utf-8") as h:
        for line in h:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 1:
                seqids.add(cols[0])
    return seqids


def filter_fastas_by_gff(gff_dir: Path, fasta_dir: Path, out_dir: Path) -> None:
    ensure_dir(out_dir)

    gffs = sorted(gff_dir.glob("*.gff"))
    fastas = list_fastas(fasta_dir)

    # Map by basename (stem)
    gff_map: Dict[str, Path] = {p.stem: p for p in gffs}
    fasta_map: Dict[str, Path] = {p.stem: p for p in fastas}

    common = sorted(set(gff_map.keys()) & set(fasta_map.keys()))
    if not common:
        print(f"WARNING: No matching basenames between GFF dir and FASTA dir:\n  {gff_dir}\n  {fasta_dir}")
        return

    for key in common:
        gff_file = gff_map[key]
        fasta_file = fasta_map[key]

        keep_ids = read_gff_seqids(gff_file)
        if not keep_ids:
            print(f"No seqids in GFF: {gff_file.name} (skipping).")
            continue

        recs = load_fasta_records(fasta_file)
        kept = [r for r in recs if fasta_id(r) in keep_ids]

        out_fa = out_dir / f"{key}_filtered.fasta"
        if kept:
            with out_fa.open("w", encoding="utf-8") as h:
                SeqIO.write(kept, h, "fasta")
            print(f"Filtered FASTA written: {out_fa}")
        else:
            print(f"No sequences kept for: {key}")


def build_blast_db(makeblastdb: str, db_fasta: Path, db_prefix: Path) -> None:
    cmd = [
        makeblastdb,
        "-in", str(db_fasta),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ]
    run_cmd(cmd, capture_stdout=False)


def blast_filtered_vs_pbs(
    filtered_dir: Path,
    pbs_fasta: Path,
    out_dir: Path,
    blastn: str,
    makeblastdb: str,
    db_prefix: Optional[Path] = None,
) -> None:
    ensure_dir(out_dir)
    if db_prefix is None:
        db_prefix = out_dir / "PBS_db"

    build_blast_db(makeblastdb, pbs_fasta, db_prefix)

    fmt = f"6 {BLAST_FIELDS}"
    for fa in list_fastas(filtered_dir):
        base = fa.stem
        out_tsv = out_dir / f"{base}_vs_PBS_db.tsv"
        cmd = [
            blastn,
            "-query", str(fa),
            "-db", str(db_prefix),
            "-out", str(out_tsv),
            "-outfmt", fmt,
        ]
        run_cmd(cmd, capture_stdout=False)
        print(f"PBS BLAST written: {out_tsv}")


# ----------------------------
# CLI
# ----------------------------

def main() -> int:
    p = argparse.ArgumentParser(
        description="Self-BLAST FASTAs in a folder, emit TSV + GFF, optionally filter FASTAs by GFF and BLAST vs PBS DB."
    )
    p.add_argument("--input-dir", "-i", required=True, type=Path, help="Directory containing input FASTA files.")
    p.add_argument("--out-dir", "-o", type=Path, default=Path("selfblast_out"), help="Output directory.")
    p.add_argument(
        "--combined-gff",
        action="store_true",
        default=False,
        help="Also write combined-coordinate GFF (<base>_selfblast_combined.gff).",
    )
    p.add_argument("--blastn", default="blastn", help="Path/name of blastn executable (default: blastn).")

    # Filtering step
    p.add_argument("--filter-fasta-dir", type=Path, default=None,
                   help="If set: directory of FASTAs to filter using produced GFFs (matched by basename).")
    p.add_argument("--filter-out-dir", type=Path, default=None,
                   help="If set: output directory for filtered FASTAs.")

    # PBS step
    p.add_argument("--pbs-fasta", type=Path, default=None, help="FASTA used to build PBS BLAST database.")
    p.add_argument("--pbs-out-dir", type=Path, default=None, help="Output directory for PBS BLAST results.")
    p.add_argument("--makeblastdb", default="makeblastdb", help="Path/name of makeblastdb executable.")

    args = p.parse_args()

    args.input_dir = args.input_dir.resolve()
    args.out_dir = args.out_dir.resolve()

    if not args.input_dir.exists():
        raise SystemExit(f"ERROR: input dir not found: {args.input_dir}")

    # Sanity: tools
    which_or_fail(args.blastn)
    if args.pbs_fasta is not None:
        which_or_fail(args.makeblastdb)

    ensure_dir(args.out_dir)

    fasta_files = list_fastas(args.input_dir)
    if not fasta_files:
        raise SystemExit(f"ERROR: No FASTA files found in {args.input_dir}")

    # Step 1: self-blast per input FASTA
    for fa in fasta_files:
        process_one_fasta(
            fasta_file=fa,
            out_dir=args.out_dir,
            blastn=args.blastn,
            write_combined_gff=args.combined_gff,
        )

    # Step 2: optional filtering (match basenames)
    filtered_dir_to_use: Optional[Path] = None
    if args.filter_fasta_dir and args.filter_out_dir:
        filter_fastas_by_gff(
            gff_dir=args.out_dir,
            fasta_dir=args.filter_fasta_dir.resolve(),
            out_dir=args.filter_out_dir.resolve(),
        )
        filtered_dir_to_use = args.filter_out_dir.resolve()

    # Step 3: optional PBS blast
    if args.pbs_fasta and args.pbs_out_dir:
        filtered_dir = filtered_dir_to_use if filtered_dir_to_use else args.out_dir
        blast_filtered_vs_pbs(
            filtered_dir=filtered_dir,
            pbs_fasta=args.pbs_fasta.resolve(),
            out_dir=args.pbs_out_dir.resolve(),
            blastn=args.blastn,
            makeblastdb=args.makeblastdb,
        )

    print("\n✅ Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())