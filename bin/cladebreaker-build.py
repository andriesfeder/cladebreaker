#!/usr/bin/env python3
"""
cladebreaker build

Builds or downloads a WhatsGNU database for a given NCBI TaxID.
If a pre-built database is available for the species, the user is prompted
to download it. Otherwise, a custom database is built from NCBI reference
genomes via the cladebreaker build Nextflow workflow.

Usage:
    cladebreaker_build.py --taxid 1280
    cladebreaker_build.py --taxid 562 --ortholog --genome_count 500
    cladebreaker_build.py --taxid 573 --outdir /data/databases -profile conda
"""

import argparse
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()


def _require_biopython():
    """Import Biopython Entrez lazily so --help always works without it installed."""
    try:
        from Bio import Entrez
        return Entrez
    except ImportError:
        print("ERROR: Biopython is required. Install with: mamba install -c bioconda biopython")
        sys.exit(1)
PROJECT_DIR = SCRIPT_DIR.parent
DATABASES_JSON = PROJECT_DIR / "assets" / "whatsgnu_databases.json"


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def get_species_name(taxid):
    """Resolve NCBI TaxID to scientific species name."""
    Entrez = _require_biopython()
    Entrez.email = "cladebreaker@pipeline.org"
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), rettype="xml")
    records = Entrez.read(handle)
    handle.close()
    if records:
        return records[0]["ScientificName"]
    return None


def count_ncbi_genomes(taxid):
    """Count assemblies available in NCBI for a given TaxID."""
    Entrez = _require_biopython()
    Entrez.email = "cladebreaker@pipeline.org"
    term = (
        f"txid{taxid}[Organism] AND latest[filter] "
        "AND all[filter] NOT anomalous[filter]"
    )
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])


# ---------------------------------------------------------------------------
# Database registry helpers
# ---------------------------------------------------------------------------

def load_registry():
    """Load the pre-built database registry from assets/whatsgnu_databases.json."""
    with open(DATABASES_JSON) as f:
        return json.load(f)["databases"]


def find_prebuilt(taxid, registry, mode):
    """Return matching pre-built database entries for taxid + mode."""
    return [db for db in registry if db["taxid"] == taxid and db["mode"] == mode]


# ---------------------------------------------------------------------------
# Download helper
# ---------------------------------------------------------------------------

def download_prebuilt(db, outdir):
    """Download a pre-built WhatsGNU database and extract if compressed."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / db["filename"]

    print(f"\nDownloading {db['filename']} to {outdir} ...")
    print("(This may take a while — database files can be large)\n")

    if subprocess.run(["which", "wget"], capture_output=True).returncode == 0:
        subprocess.run(["wget", "-c", "-O", str(outfile), db["url"]], check=True)
    else:
        urllib.request.urlretrieve(db["url"], outfile)

    # Extract compressed archives
    filename = db["filename"]
    if filename.endswith(".tar.gz"):
        print(f"Extracting {outfile} ...")
        subprocess.run(["tar", "-xzf", str(outfile), "-C", str(outdir)], check=True)
    elif filename.endswith(".zip"):
        print(f"Extracting {outfile} ...")
        subprocess.run(["unzip", "-o", str(outfile), "-d", str(outdir)], check=True)
    elif filename.endswith(".gz"):
        print(f"Decompressing {outfile} ...")
        subprocess.run(["gunzip", "-k", str(outfile)], check=True)

    print(f"\nDatabase saved to: {outdir}")


# ---------------------------------------------------------------------------
# Nextflow launch helper
# ---------------------------------------------------------------------------

def launch_build_workflow(taxid, genome_count, mode, outdir, extra_args):
    """Launch the Nextflow CLADEBREAKER_BUILD workflow for custom database building."""
    cmd = [
        "nextflow", "run", str(PROJECT_DIR / "main.nf"),
        "-entry", "CLADEBREAKER_BUILD",
        "--taxid", str(taxid),
        "--genome_count", str(genome_count),
        "--db_mode", mode,
        "--outdir", str(outdir),
    ]

    if extra_args.profile:
        cmd += ["-profile", extra_args.profile]
    if extra_args.config:
        cmd += ["-c", extra_args.config]
    if extra_args.work_dir:
        cmd += ["-w", extra_args.work_dir]

    print("\nLaunching cladebreaker build workflow:")
    print("  " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# Prompt helpers
# ---------------------------------------------------------------------------

def prompt_yes_no(question, default="yes"):
    """Prompt the user for a yes/no answer."""
    hint = "[Y/n]" if default == "yes" else "[y/N]"
    while True:
        answer = input(f"{question} {hint}: ").strip().lower()
        if answer == "" and default:
            return default == "yes"
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("  Please enter 'y' or 'n'.")


def print_banner(title):
    width = 60
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="cladebreaker build",
        description="Build or download a WhatsGNU database for a given NCBI TaxID.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  cladebreaker_build.py --taxid 1280
  cladebreaker_build.py --taxid 562 --ortholog --genome_count 500
  cladebreaker_build.py --taxid 573 --outdir /data/databases -profile conda -c cluster.config
        """,
    )

    parser.add_argument(
        "--taxid", required=True, type=int,
        help="NCBI TaxID of the target species (e.g. 1280 for S. aureus)"
    )
    parser.add_argument(
        "--genome_count", type=int, default=None,
        help="Max genomes to use for a custom build (default: all available)"
    )
    parser.add_argument(
        "--ortholog", action="store_true",
        help="Build/download an ortholog-mode database (default: basic)"
    )
    parser.add_argument(
        "--outdir", default="./cladebreaker_databases",
        help="Output directory for the database (default: ./cladebreaker_databases)"
    )
    parser.add_argument(
        "-profile", dest="profile", default=None,
        help="Nextflow profile(s) to use for a custom build (e.g. conda, docker)"
    )
    parser.add_argument(
        "-c", "--config", dest="config", default=None,
        help="Path to additional Nextflow config file for a custom build"
    )
    parser.add_argument(
        "-w", "--work-dir", dest="work_dir", default=None,
        help="Nextflow work directory for a custom build (default: ./work)"
    )

    args = parser.parse_args()
    mode = "ortholog" if args.ortholog else "basic"

    print_banner("cladebreaker build")

    # ------------------------------------------------------------------
    # Step 1: Resolve TaxID → species name
    # ------------------------------------------------------------------
    print(f"\nLooking up TaxID {args.taxid} in NCBI taxonomy ...")
    try:
        species = get_species_name(args.taxid)
    except Exception as e:
        print(f"ERROR: Could not resolve TaxID {args.taxid}: {e}")
        sys.exit(1)

    if not species:
        print(f"ERROR: TaxID {args.taxid} was not found in NCBI taxonomy.")
        sys.exit(1)

    print(f"  Species : {species}")
    print(f"  TaxID   : {args.taxid}")
    print(f"  DB mode : {mode}")

    # ------------------------------------------------------------------
    # Step 2: Check for a pre-built database
    # ------------------------------------------------------------------
    registry = load_registry()
    prebuilt = find_prebuilt(args.taxid, registry, mode)

    if prebuilt:
        db = prebuilt[0]
        size_str = f"{db['size_gb']} GB" if db["size_gb"] else "size unknown"
        print(f"\nPre-built {mode} database found for {species}:")
        print(f"  Version  : {db['version']}")
        print(f"  Genomes  : {db['genome_count']:,}")
        print(f"  Size     : {size_str}")

        use_prebuilt = prompt_yes_no(
            f"\nWould you like to download this pre-built database?"
        )

        if use_prebuilt:
            print(f"\nDownload location : {args.outdir}")
            if prompt_yes_no("Proceed with download?"):
                download_prebuilt(db, args.outdir)
                print("\nDone! Point cladebreaker to this database with --db")
            else:
                print("Download cancelled.")
            return

        print("\nSkipping pre-built database. Proceeding to custom build.")

    # ------------------------------------------------------------------
    # Step 3: Count available NCBI genomes
    # ------------------------------------------------------------------
    print(f"\nQuerying NCBI Assembly for {species} (TaxID: {args.taxid}) ...")
    try:
        available = count_ncbi_genomes(args.taxid)
    except Exception as e:
        print(f"ERROR: Could not query NCBI Assembly database: {e}")
        sys.exit(1)

    if available == 0:
        print(f"ERROR: No assemblies found in NCBI for TaxID {args.taxid}.")
        sys.exit(1)

    genome_count = args.genome_count or available
    if genome_count > available:
        print(
            f"  Warning: Requested {genome_count:,} genomes but only "
            f"{available:,} are available. Using {available:,}."
        )
        genome_count = available

    # ------------------------------------------------------------------
    # Step 4: Confirmation prompt
    # ------------------------------------------------------------------
    print_banner("Build Configuration")
    print(f"  Species          : {species}")
    print(f"  TaxID            : {args.taxid}")
    print(f"  Genomes to use   : {genome_count:,}  (of {available:,} available in NCBI)")
    print(f"  Database mode    : {mode}")
    print(f"  Output directory : {args.outdir}")
    if args.profile:
        print(f"  Nextflow profile : {args.profile}")
    if args.config:
        print(f"  Nextflow config  : {args.config}")
    if args.work_dir:
        print(f"  Nextflow work dir: {args.work_dir}")
    print("=" * 60)

    if not prompt_yes_no("\nProceed with build?"):
        print("Build cancelled.")
        sys.exit(0)

    # ------------------------------------------------------------------
    # Step 5: Launch Nextflow build workflow
    # ------------------------------------------------------------------
    launch_build_workflow(args.taxid, genome_count, mode, args.outdir, args)


if __name__ == "__main__":
    main()
