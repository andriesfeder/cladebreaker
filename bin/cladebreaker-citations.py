#!/usr/bin/env python3
"""
cladebreaker citations

Prints citations for cladebreaker and all tools it uses.
"""

import argparse
import textwrap


CITATIONS = [
    {
        "tool": "cladebreaker",
        "citation": (
            "Feder AF, Planet PJ (2022). cladebreaker: a Nextflow pipeline for "
            "whole-genome sequencing-based clade assignment of bacterial isolates. "
            "https://github.com/andriesfeder/cladebreaker"
        ),
    },
    {
        "tool": "WhatsGNU",
        "citation": (
            "Moustafa AM, Planet PJ (2020). WhatsGNU: a tool for identifying "
            "proteomic novelty. Genome Biology 21:58. "
            "https://doi.org/10.1186/s13059-020-01965-w"
        ),
    },
    {
        "tool": "Nextflow",
        "citation": (
            "Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C "
            "(2017). Nextflow enables reproducible computational workflows. "
            "Nature Biotechnology 35:316-319. "
            "https://doi.org/10.1038/nbt.3820"
        ),
    },
    {
        "tool": "FastQC",
        "citation": (
            "Andrews S (2010). FastQC: a quality control tool for high throughput "
            "sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
        ),
    },
    {
        "tool": "Shovill",
        "citation": (
            "Seemann T (2019). Shovill: faster SPAdes assembly of Illumina reads. "
            "https://github.com/tseemann/shovill"
        ),
    },
    {
        "tool": "Prokka",
        "citation": (
            "Seemann T (2014). Prokka: rapid prokaryotic genome annotation. "
            "Bioinformatics 30:2068-2069. "
            "https://doi.org/10.1093/bioinformatics/btu153"
        ),
    },
    {
        "tool": "Roary",
        "citation": (
            "Page AJ, Cummins CA, Hunt M, Wong VK, Reuter S, Holden MTG, Fookes M, "
            "Falush D, Keane JA, Parkhill J (2015). Roary: rapid large-scale prokaryote "
            "pan genome analysis. Bioinformatics 31:3691-3693. "
            "https://doi.org/10.1093/bioinformatics/btv421"
        ),
    },
    {
        "tool": "PIRATE",
        "citation": (
            "Bayliss SC, Thorpe HA, Coyle NM, Sheppard SK, Feil EJ (2019). PIRATE: "
            "a fast and scalable pangenomics toolbox for clustering diverged orthologues "
            "in bacteria. GigaScience 8:giz119. "
            "https://doi.org/10.1093/gigascience/giz119"
        ),
    },
    {
        "tool": "RAxML-NG",
        "citation": (
            "Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A (2019). RAxML-NG: "
            "a fast, scalable, and user-friendly tool for maximum likelihood phylogenetic "
            "inference. Bioinformatics 35:4453-4455. "
            "https://doi.org/10.1093/bioinformatics/btz305"
        ),
    },
    {
        "tool": "Snippy",
        "citation": (
            "Seemann T (2015). Snippy: fast bacterial variant calling from NGS reads. "
            "https://github.com/tseemann/snippy"
        ),
    },
    {
        "tool": "assembly-scan",
        "citation": (
            "Petit RA III (2019). assembly-scan: generate basic stats for an assembly. "
            "https://github.com/rpetit3/assembly-scan"
        ),
    },
    {
        "tool": "ncbi-genome-download",
        "citation": (
            "Blin K et al. ncbi-genome-download: scripts to download genomes from "
            "NCBI's FTP servers. https://github.com/kblin/ncbi-genome-download"
        ),
    },
    {
        "tool": "MultiQC",
        "citation": (
            "Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize "
            "analysis results for multiple tools and samples in a single report. "
            "Bioinformatics 32:3047-3048. "
            "https://doi.org/10.1093/bioinformatics/btw354"
        ),
    },
]


def print_citations(tool_filter=None):
    width = 80
    sep = "-" * width

    tools = CITATIONS
    if tool_filter:
        tools = [c for c in CITATIONS if tool_filter.lower() in c["tool"].lower()]
        if not tools:
            print(f"No citations found matching '{tool_filter}'.")
            return

    print()
    for entry in tools:
        print(f"  [{entry['tool']}]")
        wrapped = textwrap.fill(entry["citation"], width=width - 4, initial_indent="    ", subsequent_indent="    ")
        print(wrapped)
        print(sep)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="cladebreaker citations",
        description="Print citations for cladebreaker and all tools it uses.",
    )
    parser.add_argument(
        "--tool", metavar="NAME", default=None,
        help="Filter citations to a specific tool by name.",
    )
    args = parser.parse_args()
    print_citations(args.tool)


if __name__ == "__main__":
    main()
