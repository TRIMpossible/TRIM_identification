#!/usr/bin/env python3
"""
Runs scripts:

  00_filter_BLAST_by_query.py
  01_convert_blast_2_gff.py
  02_merge_gff_entries.py
  03_fas_from_gff.py
  04_selfblast_pipeline.py

What is needed?
-----------
- You have the 5 scripts in the SAME folder as this master script:
    quantification_output_precessing.py
    00_filter_BLAST_by_query.py
    01_convert_blast_2_gff.py
    02_merge_gff_entries.py
    03_fas_from_gff.py
    selfblast_pipeline.py

- BLAST+ tools in PATH: blastn (and optionally makeblastdb if you use PBS later)
- Step 03 dependencies installed if you run it:
    pip install biopython bcbio-gff
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import shutil


def die(msg: str, code: int = 2) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def run(cmd: list[str]) -> None:
    print("\n>> " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def which_or_die(prog: str) -> None:
    if not shutil.which(prog):
        die(f"Required executable not found in PATH: {prog}")


def must_exist(p: Path, what: str) -> None:
    if not p.exists():
        die(f"{what} not found: {p}")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Master pipeline: BLAST TSV -> per-query -> GFF -> merged GFF -> FASTA (split/non-split) -> self-BLAST."
    )

    # Required inputs
    ap.add_argument("--blast", "-b", required=True, type=Path, help="Input BLAST tabular file (outfmt 6/7-like).")
    ap.add_argument("--genome", "-g", required=True, type=Path, help="Genome/reference FASTA file.")
    ap.add_argument("--workdir", "-w", default=Path("master_out"), type=Path, help="Output base directory.")

    # Step 00 options
    ap.add_argument("--query-col", type=int, default=0, help="0-based query column index for step 00 (default: 0).")
    ap.add_argument("--queries", type=Path, default=None, help="Optional file of query IDs to keep (step 00).")
    ap.add_argument("--split-suffix", default="_filtered_output.txt", help="Suffix used by step 00 outputs.")

    # Step 01 options
    ap.add_argument("--blast-pattern", default="*.txt", help="Pattern for step 01 input files in step00 dir (default: *.txt).")
    ap.add_argument("--family-tag-from", choices=["prefix", "filename"], default="prefix",
                    help="Step 01 family tag derivation (default: prefix).")

    # Step 02 options
    ap.add_argument("--merge-distance", type=int, default=5000, help="Step 02 merge distance (bp).")
    ap.add_argument("--identity-threshold", type=float, default=75.0, help="Step 02 identity threshold (%).")
    ap.add_argument("--split-threshold", type=int, default=10000, help="Split threshold (bp) used by step 02 and 03.")

    # Step 03 options
    ap.add_argument("--flanks", default="2000", help="Comma-separated flank sizes for step 03 (default: 2000).")
    ap.add_argument("--gff-pattern", default="*_merged.gff", help="Pattern for step 03 input GFFs (default: *_merged.gff).")

    # Selfblast options
    ap.add_argument("--run-selfblast", action="store_true", help="Run selfblast_pipeline.py on the split FASTAs.")
    ap.add_argument("--selfblast-combined-gff", action="store_true",
                    help="Pass --combined-gff to selfblast_pipeline.py (writes combined coordinate GFF).")

    # IMPORTANT: grouping in selfblast script
    ap.add_argument(
        "--selfblast-input",
        choices=["split", "non_split"],
        default="split",
        help="Which FASTAs from step 03 to self-BLAST: split or non_split (default: split).",
    )

    args = ap.parse_args()

    # Resolve & sanity check inputs
    args.blast = args.blast.resolve()
    args.genome = args.genome.resolve()
    args.workdir = args.workdir.resolve()

    must_exist(args.blast, "BLAST input file")
    must_exist(args.genome, "Genome FASTA")

    # Scripts (must be in same folder as master)
    here = Path(__file__).resolve().parent
    py = sys.executable

    s00 = here / "00_filter_BLAST_by_query.py"
    s01 = here / "01_convert_blast_2_gff.py"
    s02 = here / "02_merge_gff_entries.py"
    s03 = here / "03_fas_from_gff.py"
    sSB = here / "selfblast_pipeline.py"

    for s, name in [(s00, "step00"), (s01, "step01"), (s02, "step02"), (s03, "step03")]:
        must_exist(s, f"Script {name}")

    if args.run_selfblast:
        must_exist(sSB, "Script selfblast_pipeline.py")
        which_or_die("blastn")

    # Create output structure
    step00_dir = args.workdir / "00_filtered"
    step01_dir = args.workdir / "01_gff"
    step02_dir = args.workdir / "02_merged"
    logs_dir = args.workdir / "02_logs"
    step03_dir = args.workdir / "03_fasta"

    for d in [step00_dir, step01_dir, step02_dir, logs_dir, step03_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # -----------------------
    # Step 00
    # -----------------------
    cmd00 = [
        py, str(s00),
        "-i", str(args.blast),
        "-o", str(step00_dir),
        "--query-col", str(args.query_col),
        "--suffix", args.split_suffix,
    ]
    if args.queries:
        args.queries = args.queries.resolve()
        must_exist(args.queries, "Queries file")
        cmd00 += ["--queries", str(args.queries)]
    run(cmd00)

    # -----------------------
    # Step 01
    # -----------------------
    cmd01 = [
        py, str(s01),
        "-i", str(step00_dir),
        "-o", str(step01_dir),
        "--pattern", args.blast_pattern,
        "--family-tag-from", args.family_tag_from,
    ]
    run(cmd01)

    # -----------------------
    # Step 02
    # -----------------------
    cmd02 = [
        py, str(s02),
        "-i", str(step01_dir),
        "-o", str(step02_dir),
        "--log-dir", str(logs_dir),
        "--merge-distance", str(args.merge_distance),
        "--identity-threshold", str(args.identity_threshold),
        "--split-threshold", str(args.split_threshold),
    ]
    run(cmd02)

    # -----------------------
    # Step 03
    # -----------------------
    cmd03 = [
        py, str(s03),
        "-i", str(step02_dir),
        "-g", str(args.genome),
        "-o", str(step03_dir),
        "--pattern", args.gff_pattern,
        "--flanks", args.flanks,
        "--split-threshold", str(args.split_threshold),
    ]
    run(cmd03)

    # -----------------------
    # Selfblast (optional)
    # -----------------------
    if args.run_selfblast:
        # Step03 writes into subfolders: flank<flanksize>/
        # We pick the FIRST flank size as default target for selfblast (common use-case).
        flank_sizes = [x.strip() for x in str(args.flanks).split(",") if x.strip()]
        if not flank_sizes:
            die("No valid --flanks provided (e.g. 500,2000).")
        flank0 = flank_sizes[0]
        flank_dir = step03_dir / f"flank{flank0}"
        if not flank_dir.exists():
            die(f"Expected flank directory not found: {flank_dir}")

        # selfblast needs input FASTAs in a directory. We create a directory with only the chosen FASTAs.
        sb_in = args.workdir / f"04_selfblast_input_flank{flank0}_{args.selfblast_input}"
        sb_out = args.workdir / f"04_selfblast_out_flank{flank0}_{args.selfblast_input}"
        sb_in.mkdir(parents=True, exist_ok=True)
        sb_out.mkdir(parents=True, exist_ok=True)

        pattern = "*_split.fasta" if args.selfblast_input == "split" else "*_non_split.fasta"
        fastas = sorted(flank_dir.glob(pattern))
        if not fastas:
            die(f"No FASTA files matched {pattern} in {flank_dir}")

        # Copy (not symlink) for Windows-friendliness
        for fa in fastas:
            shutil.copy2(fa, sb_in / fa.name)

        cmdSB = [
            py, str(sSB),
            "--input-dir", str(sb_in),
            "--out-dir", str(sb_out),
        ]
        if args.selfblast_combined_gff:
            cmdSB.append("--combined-gff")
        run(cmdSB)

        print(f"\nSelfblast finished. Input copied to: {sb_in}")
        print(f"Selfblast outputs in: {sb_out}")

    print(f"\n✅ Master pipeline finished. Outputs under: {args.workdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())