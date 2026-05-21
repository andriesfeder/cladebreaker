#!/usr/bin/env python3
"""
cladebreaker prepare

Generates a cladebreaker input samplesheet (CSV) from a directory of FASTQ
files and/or genome assemblies. Paired-end reads are matched automatically
using common R1/R2 filename conventions.

Output format:
    sample,type,file_1,file_2

Usage:
    cladebreaker prepare --fastqs /path/to/reads/ --output samplesheet.csv
    cladebreaker prepare --assemblies /path/to/assemblies/ --output samplesheet.csv
    cladebreaker prepare --fastqs /path/to/reads/ --assemblies /path/to/assemblies/ --output samplesheet.csv
"""

import argparse
import os
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
ASSEMBLY_EXTENSIONS = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")

# Patterns used to strip R1/R2 suffixes and detect read direction.
# Listed from most specific to least specific.
R1_PATTERNS = [
    r"_R1_\d+",   # _R1_001
    r"_R1\b",     # _R1
    r"_1\b",      # _1
]
R2_PATTERNS = [
    r"_R2_\d+",   # _R2_001
    r"_R2\b",     # _R2
    r"_2\b",      # _2
]


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------

def find_files(directory, extensions, recursive=False):
    """Return all files in directory matching any of the given extensions."""
    directory = Path(directory)
    pattern = "**/*" if recursive else "*"
    files = []
    for ext in extensions:
        files.extend(directory.glob(f"{pattern}{ext}"))
    return sorted(set(files))


def strip_read_suffix(name, patterns):
    """Remove a read-direction suffix from a filename stem and return (stem, matched)."""
    for pat in patterns:
        m = re.search(pat, name)
        if m:
            stem = name[: m.start()]
            return stem, True
    return name, False


def pair_fastqs(files):
    """
    Given a list of FASTQ paths, return a list of (sample, r1, r2) tuples.
    r2 is None for single-end files.
    """
    r1_map = {}   # sample_stem -> Path
    r2_map = {}   # sample_stem -> Path
    unpaired = [] # files that don't match any R1/R2 pattern

    for f in files:
        stem = f.name
        for ext in FASTQ_EXTENSIONS:
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break

        base, is_r1 = strip_read_suffix(stem, R1_PATTERNS)
        if is_r1:
            r1_map[base] = f
            continue

        base, is_r2 = strip_read_suffix(stem, R2_PATTERNS)
        if is_r2:
            r2_map[base] = f
            continue

        unpaired.append(f)

    samples = []

    # Paired-end
    for base in sorted(r1_map):
        r1 = r1_map[base]
        r2 = r2_map.get(base)
        if r2:
            samples.append((base, "reads", str(r1.resolve()), str(r2.resolve())))
        else:
            # R1 with no matching R2 — treat as single-end
            samples.append((base, "reads", str(r1.resolve()), ""))

    # R2 with no matching R1 (unusual, include as single-end with warning)
    for base in sorted(r2_map):
        if base not in r1_map:
            print(
                f"  WARNING: Found R2 file with no matching R1 for sample '{base}' — skipping.",
                file=sys.stderr,
            )

    # Files that matched no R1/R2 pattern — treat as single-end
    for f in unpaired:
        stem = f.name
        for ext in FASTQ_EXTENSIONS:
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        samples.append((stem, "reads", str(f.resolve()), ""))

    return samples


def find_assemblies(files):
    """Return a list of (sample, type, file_1, file_2) for assembly files."""
    samples = []
    for f in files:
        stem = f.name
        for ext in ASSEMBLY_EXTENSIONS:
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        samples.append((stem, "assembly", str(f.resolve()), ""))
    return samples


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_samplesheet(samples, output):
    """Write samples to a CSV samplesheet."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    deduped = []
    for row in samples:
        key = row[0]
        if key in seen:
            print(
                f"  WARNING: Duplicate sample name '{key}' detected — keeping first occurrence.",
                file=sys.stderr,
            )
            continue
        seen.add(key)
        deduped.append(row)

    with open(output, "w") as fh:
        fh.write("sample,type,file_1,file_2\n")
        for sample, stype, f1, f2 in sorted(deduped, key=lambda x: x[0]):
            fh.write(f"{sample},{stype},{f1},{f2}\n")

    return len(deduped)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="cladebreaker prepare",
        description=(
            "Generate a cladebreaker input samplesheet from a directory of "
            "FASTQ files and/or genome assemblies."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cladebreaker prepare --fastqs /data/reads/ --output samplesheet.csv
  cladebreaker prepare --assemblies /data/assemblies/ --output samplesheet.csv
  cladebreaker prepare --fastqs /data/reads/ --assemblies /data/assemblies/ --output samplesheet.csv
  cladebreaker prepare --fastqs /data/reads/ --recursive --output samplesheet.csv
        """,
    )

    parser.add_argument(
        "--fastqs", metavar="DIR", default=None,
        help="Directory of FASTQ files (.fastq.gz, .fq.gz). Paired-end files are matched automatically."
    )
    parser.add_argument(
        "--assemblies", metavar="DIR", default=None,
        help="Directory of genome assembly files (.fna, .fa, .fasta)."
    )
    parser.add_argument(
        "--output", "-o", metavar="FILE", required=True,
        help="Path to write the output samplesheet CSV."
    )
    parser.add_argument(
        "--recursive", "-r", action="store_true",
        help="Search for files recursively within the input directory/directories."
    )
    parser.add_argument(
        "--prefix", metavar="STR", default=None,
        help="Only include files whose names begin with this prefix."
    )

    args = parser.parse_args()

    if not args.fastqs and not args.assemblies:
        parser.error("At least one of --fastqs or --assemblies must be provided.")

    all_samples = []

    # Process FASTQ directory
    if args.fastqs:
        fastq_dir = Path(args.fastqs)
        if not fastq_dir.is_dir():
            print(f"ERROR: --fastqs directory does not exist: {fastq_dir}", file=sys.stderr)
            sys.exit(1)

        fastq_files = find_files(fastq_dir, FASTQ_EXTENSIONS, args.recursive)
        if args.prefix:
            fastq_files = [f for f in fastq_files if f.name.startswith(args.prefix)]

        if not fastq_files:
            print(f"WARNING: No FASTQ files found in {fastq_dir}", file=sys.stderr)
        else:
            paired = pair_fastqs(fastq_files)
            all_samples.extend(paired)
            pe = sum(1 for s in paired if s[3])
            se = sum(1 for s in paired if not s[3])
            print(f"Found {len(fastq_files)} FASTQ files → {pe} paired-end, {se} single-end samples.")

    # Process assembly directory
    if args.assemblies:
        asm_dir = Path(args.assemblies)
        if not asm_dir.is_dir():
            print(f"ERROR: --assemblies directory does not exist: {asm_dir}", file=sys.stderr)
            sys.exit(1)

        asm_files = find_files(asm_dir, ASSEMBLY_EXTENSIONS, args.recursive)
        if args.prefix:
            asm_files = [f for f in asm_files if f.name.startswith(args.prefix)]

        if not asm_files:
            print(f"WARNING: No assembly files found in {asm_dir}", file=sys.stderr)
        else:
            asms = find_assemblies(asm_files)
            all_samples.extend(asms)
            print(f"Found {len(asm_files)} assembly files → {len(asms)} samples.")

    if not all_samples:
        print("ERROR: No samples found. Check your input directories and extensions.", file=sys.stderr)
        sys.exit(1)

    count = write_samplesheet(all_samples, args.output)
    print(f"Samplesheet written to {args.output} ({count} samples).")
    print(f"\nRun cladebreaker with:")
    print(f"  cladebreaker --input {args.output} --db /path/to/db.pickle [options]")


if __name__ == "__main__":
    main()
