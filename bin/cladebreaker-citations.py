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
        "tool": "Panaroo",
        "citation": (
            "Tonkin-Hill G, MacAlasdair N, Gillies C, Weimann A, Rocha EPC, Bentley SD, "
            "Parkhill J, Page AJ, Getz A, Wordsworth S, Coll F (2020). Producing polished "
            "prokaryotic pangenomes with the Panaroo pipeline. Genome Biology 21:180. "
            "https://doi.org/10.1186/s13059-020-02090-4"
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
        "tool": "Bakta",
        "citation": (
            "Schwengers O, Jelonek L, Dieckmann MA, Beyvers S, Blom J, Goesmann A (2021). "
            "Bakta: rapid and standardized annotation of bacterial genomes via "
            "alignment-free sequence identification. Microbial Genomics 7:000685. "
            "https://doi.org/10.1099/mgen.0.000685"
        ),
    },
    {
        "tool": "genealogical sorting index",
        "citation": (
            "Cummings MP, Neel MC, Shaw KL (2008). A genealogical approach to "
            "quantifying lineage divergence. Evolution 62:2411-2422. "
            "https://doi.org/10.1111/j.1558-5646.2008.00442.x"
        ),
    },
    {
        "tool": "Rosenberg's tests for chance monophyly",
        "citation": (
            "Rosenberg NA (2007). Statistical tests for taxonomic distinctiveness "
            "from observations of monophyly. Evolution 61:317-323. "
            "https://doi.org/10.1111/j.1558-5646.2007.00023.x"
        ),
    },
    {
        "tool": "joint monophyly of k groups",
        "citation": (
            "Zhu S, Degnan JH, Steel M (2011). Clades, clans and reciprocal monophyly "
            "under neutral evolutionary models. Theoretical Population Biology 80:118-124. "
            "https://doi.org/10.1016/j.tpb.2011.05.002"
        ),
    },
    {
        "tool": "Slatkin-Maddison migration test",
        "citation": (
            "Slatkin M, Maddison WP (1989). A cladistic measure of gene flow inferred "
            "from the phylogenies of alleles. Genetics 123:603-613. "
            "https://doi.org/10.1093/genetics/123.3.603"
        ),
    },
    {
        "tool": "snp-dists",
        "citation": (
            "Seemann T. snp-dists: pairwise SNP distance matrix from a FASTA "
            "sequence alignment. https://github.com/tseemann/snp-dists"
        ),
    },
    {
        "tool": "DendroPy",
        "citation": (
            "Sukumaran J, Holder MT (2010). DendroPy: a Python library for "
            "phylogenetic computing. Bioinformatics 26:1569-1571. "
            "https://doi.org/10.1093/bioinformatics/btq228"
        ),
    },
    {
        "tool": "Biopython",
        "citation": (
            "Cock PJ, Antao T, Chang JT, et al. (2009). Biopython: freely available "
            "Python tools for computational molecular biology and bioinformatics. "
            "Bioinformatics 25:1422-1423. https://doi.org/10.1093/bioinformatics/btp163"
        ),
    },
    {
        "tool": "Matplotlib",
        "citation": (
            "Hunter JD (2007). Matplotlib: A 2D graphics environment. "
            "Computing in Science & Engineering 9:90-95. "
            "https://doi.org/10.1109/MCSE.2007.55"
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
