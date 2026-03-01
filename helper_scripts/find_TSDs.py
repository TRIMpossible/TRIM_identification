"""
find TSD sequences based on full length "element" annotation"
"""

from Bio import SeqIO
import sys
import os
import re
from collections import defaultdict

def parse_gff(gff_file):
    element_coords = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            if feature == "element":
                if seqid not in element_coords:
                    element_coords[seqid] = []
                element_coords[seqid].append((int(start), int(end)))
    print(f"[DEBUG] Parsed {sum(len(v) for v in element_coords.values())} 'element' features from GFF.")
    return element_coords

def get_family_name(seqid):
    match = re.match(r'(TRIM-\d+)', seqid)
    return match.group(1) if match else "Unknown"

def get_tsd_candidates(seq, start, end):
    flank5_left = seq[start-6:start-1] if start >= 6 else ""
    flank5_right = seq[end:end+5] if end+5 <= len(seq) else ""
    return flank5_left.upper(), flank5_right.upper()

def matches_with_one_mismatch(seq1, seq2):
    if len(seq1) != 5 or len(seq2) != 5:
        return False
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches <= 1

def find_tsds(gff_file, fasta_file, output_prefix="TSD", master_output=True):
    print(f"[DEBUG] Loading GFF: {gff_file}")
    element_coords = parse_gff(gff_file)

    print(f"[DEBUG] Loading FASTA: {fasta_file}")
    records = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"[DEBUG] Loaded {len(records)} FASTA records.")
    
    family_tsd_annotations = defaultdict(list)
    all_annotations = []

    for seqid in element_coords:
        matching_record = next((r for r in records if r.id == seqid), None)
        if not matching_record:
            print(f"[WARNING] Sequence '{seqid}' from GFF not found in FASTA.")
            continue

        family = get_family_name(seqid)
        sequence = str(matching_record.seq)

        for start, end in element_coords[seqid]:
            flank5_left, flank5_right = get_tsd_candidates(sequence, start, end)
            print(f"[DEBUG] Checking {seqid} ({family}): {start}-{end} | Left: '{flank5_left}' | Right: '{flank5_right}'")

            if matches_with_one_mismatch(flank5_left, flank5_right):
                print(f"[INFO] TSD found in {seqid}: {flank5_left} ~ {flank5_right}")
                # Left TSD annotation
                left_entry = (
                    seqid,
                    "TSD_finder",
                    "TSD_left",
                    start - 5,
                    start - 1,
                    ".",
                    "+",
                    ".",
                    f"Note=TSD left with <=1 mismatch;left={flank5_left};right={flank5_right}"
                )
                # Right TSD annotation
                right_entry = (
                    seqid,
                    "TSD_finder",
                    "TSD_right",
                    end + 1,
                    end + 5,
                    ".",
                    "+",
                    ".",
                    f"Note=TSD right with <=1 mismatch;left={flank5_left};right={flank5_right}"
                )

                family_tsd_annotations[family].extend([left_entry, right_entry])
                all_annotations.extend([left_entry, right_entry])

    # Write one GFF file per family
    for family, entries in family_tsd_annotations.items():
        output_file = f"{output_prefix}_{family}.gff"
        print(f"[DEBUG] Writing {len(entries)} entries to {output_file}")
        with open(output_file, "w") as out:
            out.write("##gff-version 3\n")
            for entry in entries:
                out.write("\t".join(map(str, entry)) + "\n")

    # Optionally write a master GFF
    if master_output and all_annotations:
        master_file = f"{output_prefix}_master.gff"
        print(f"[DEBUG] Writing master GFF with {len(all_annotations)} entries to {master_file}")
        with open(master_file, "w") as out:
            out.write("##gff-version 3\n")
            for entry in all_annotations:
                out.write("\t".join(map(str, entry)) + "\n")

    print(f"[DONE] Wrote {len(family_tsd_annotations)} family GFF files and master GFF.")

# Add CLI support
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python find_TSDs.py <input.gff> <input.fasta> [output_prefix]")
        sys.exit(1)
    
    gff_input = sys.argv[1]
    fasta_input = sys.argv[2]
    output_prefix = sys.argv[3] if len(sys.argv) > 3 else "TSD"

    find_tsds(gff_input, fasta_input, output_prefix)
