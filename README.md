# Cladebreaker

[![GitHub release](https://img.shields.io/github/release/andriesfeder/cladebreaker.svg)](https://github.com/andriesfeder/cladebreaker/releases)
[![Conda](https://img.shields.io/conda/vn/bioconda/cladebreaker.svg)](https://anaconda.org/bioconda/cladebreaker)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Cladebreaker** is a Nextflow pipeline for whole-genome sequencing (WGS)-based phylogenetic analysis of bacterial isolates. Starting from raw Illumina reads or assembled genomes, Cladebreaker identifies the most closely related publicly available genomes using [WhatsGNU](https://github.com/ahmedmagds/WhatsGNU), downloads those reference genomes, and builds a phylogenetic context — either through pangenome alignment ([Panaroo](https://github.com/gtonkinhill/panaroo), [PIRATE](https://github.com/SionBayliss/PIRATE), or [Roary](https://github.com/sanger-pathogens/Roary)) or reference-based SNP alignment ([Snippy](https://github.com/tseemann/snippy)) — with an optional maximum-likelihood tree via [RAxML-NG](https://github.com/amkozlov/raxml-ng).

Cladebreaker was designed to be practical: it works on any bacterial species, runs on laptops and HPC clusters alike, and requires nothing beyond a mamba install to get started.

---

## Features

- **Species-agnostic**: works for any bacterial species with NCBI assemblies available
- **Flexible input**: accepts paired-end reads, single-end reads, pre-assembled genomes, or any combination
- **WhatsGNU-powered genome selection**: identifies the closest relatives among thousands of public genomes using proteomic rarity scoring — no arbitrary species list required
- **Three pangenome tools**: Panaroo (default), PIRATE, or Roary — choose the tool that fits your analysis
- **Reference-based alignment**: use `--ref` to switch from pangenome mode to Snippy SNP alignment for any target reference
- **Integrated de novo assembly**: Shovill assembles raw reads in-pipeline, so you never have to run an assembler separately
- **Automated database management**: `cladebreaker build` downloads pre-built WhatsGNU databases for common species or constructs a custom database from NCBI for any taxon
- **Pre-built databases available**: ready-to-use databases for *S. aureus*, *S. epidermidis*, *K. pneumoniae*, *P. aeruginosa*, *E. coli*, *S. enterica*, *M. tuberculosis*, and *C. difficile*
- **Reproducible**: built on [Nextflow](https://www.nextflow.io/) with full support for Conda, Docker, Singularity, Podman, and other container backends
- **MultiQC report**: a single HTML quality-control report covering read quality, assembly statistics, and software versions

---

## Quick Start

```sh
# 1. Install
mamba create -n cladebreaker -c bioconda cladebreaker
mamba activate cladebreaker

# 2. Download a database (e.g., S. aureus ortholog)
cladebreaker build --taxid 1280

# 3. Generate a samplesheet from your reads
cladebreaker prepare --fastqs /path/to/reads/ --output samplesheet.csv

# 4. Run the pipeline
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/database.pickle \
  --o \
  -profile conda
```

That's it. Cladebreaker handles assembly, annotation, closest-genome selection, genome download, alignment, and QC reporting automatically.

---

## Installation

Cladebreaker is distributed through [Bioconda](https://bioconda.github.io/) and is best installed with [Mamba](https://github.com/mamba-org/mamba) for fast dependency resolution:

```sh
mamba create -n cladebreaker -c bioconda cladebreaker
mamba activate cladebreaker
```

Conda also works if Mamba is not available:

```sh
conda create -n cladebreaker -c bioconda cladebreaker
conda activate cladebreaker
```

Verify the installation:

```sh
cladebreaker version
# cladebreaker v0.3.1
```

> **Nextflow requirement**: Cladebreaker requires Nextflow 26.x or later. If Nextflow is not already installed in your environment, it will be included as a dependency through Bioconda.

---

## Overview

Cladebreaker is invoked through a single entry point with several subcommands:

```
cladebreaker [options]            Run the main analysis pipeline
cladebreaker build [options]      Download or build a WhatsGNU database
cladebreaker prepare [options]    Generate an input samplesheet
cladebreaker citations            Print citations for all tools
cladebreaker version              Print the cladebreaker version
```

The typical workflow follows three steps:

```
Step 1 ── cladebreaker build    →  obtain a WhatsGNU database for your species
Step 2 ── cladebreaker prepare  →  generate a samplesheet from your input files
Step 3 ── cladebreaker          →  run the full analysis pipeline
```

---

## Step 1 — Build or Download a Database (`cladebreaker build`)

Cladebreaker requires a WhatsGNU database to identify which publicly available genomes are most closely related to your query isolates. The `build` subcommand handles this automatically: it checks whether a pre-built database is available for your species, and if so prompts you to download it. If no pre-built database exists, it constructs one from scratch using NCBI reference genomes.

```sh
cladebreaker build --taxid 1280
```

The command will query NCBI taxonomy to resolve the species name, check the pre-built database registry, display a summary, and prompt you before downloading or building anything.

### Pre-built databases

Pre-built databases are available for the following species and can be downloaded interactively:

| Species | TaxID | Mode | Genomes | Size | Version |
|---|---|---|---:|---:|---|
| *Staphylococcus aureus* | 1280 | ortholog | 68,299 | 14 GB | April 2024 |
| *Clostridioides difficile* | 1496 | ortholog | 14,186 | 3.8 GB | July 2024 |
| *Klebsiella pneumoniae* | 573 | ortholog | 8,752 | — | April 2023 |
| *Pseudomonas aeruginosa* | 287 | ortholog | 4,712 | — | July 2019 |
| *Mycobacterium tuberculosis* | 1773 | ortholog | 6,563 | — | July 2019 |
| *Klebsiella pneumoniae* | 573 | basic | 75,246 | 37 GB | June 2024 |
| *Pseudomonas aeruginosa* | 287 | basic | 31,832 | 19 GB | June 2024 |
| *Escherichia coli* | 562 | basic | 211,942 | 90 GB | March 2024 |
| *Salmonella enterica* | 28901 | basic | 216,642 | — | August 2019 |
| *Staphylococcus epidermidis* | 1282 | basic | 4,981 | 0.8 GB | January 2025 |

For any species not listed, `cladebreaker build` will automatically download reference genomes from NCBI Assembly and construct a custom database.

### Database modes

WhatsGNU databases come in two modes:

- **basic**: compresses all proteins across all genomes. Suitable for any species. Recommended as the default starting point.
- **ortholog**: groups proteins by ortholog clusters before compression, providing additional rarity scoring resolution. Requires a well-sampled reference set. Best suited for species with large, curated databases.

The mode used to build the database must match the `--o` (ortholog) or `--b` (basic) flag when running the main pipeline.

### `cladebreaker build` options

| Option | Description |
|---|---|
| `--taxid INT` | NCBI TaxID of the target species (required) |
| `--ortholog` | Download/build an ortholog-mode database (default: basic) |
| `--genome_count INT` | Maximum genomes for a custom build (default: all available in NCBI) |
| `--outdir PATH` | Where to save the database (default: `./cladebreaker_databases`) |
| `-profile STR` | Nextflow profile for a custom build (e.g. `conda`, `docker`) |
| `-c PATH` | Additional Nextflow config file for a custom build |
| `-w PATH` | Nextflow work directory for a custom build |

### Examples

```sh
# Download the pre-built basic database for K. pneumoniae (75,246 genomes)
cladebreaker build --taxid 573

# Download the pre-built ortholog database for S. aureus (68,299 genomes)
cladebreaker build --taxid 1280 --ortholog

# Build a custom database for A. baumannii using 200 NCBI reference genomes
cladebreaker build --taxid 470 --genome_count 200 --outdir /data/databases -profile conda

# Build a custom ortholog database for any species, using all available genomes
cladebreaker build --taxid 1313 --ortholog --outdir /data/databases -profile docker
```

### How custom builds work

When no pre-built database is available (or you choose to build your own), the `build` subcommand launches a Nextflow workflow that:

1. Fetches genome accession lists from NCBI Assembly for the given TaxID
2. Downloads all genome FASTA files using ncbi-genome-download
3. Annotates each genome with Prokka to generate protein FAA files
4. Compresses all FAA files into a WhatsGNU database using the specified mode

---

## Step 2 — Prepare a Samplesheet (`cladebreaker prepare`)

Cladebreaker requires a CSV samplesheet describing your input samples. The `prepare` subcommand generates one automatically from a directory of FASTQ files and/or assembled genomes — no manual editing required.

```sh
cladebreaker prepare --fastqs /path/to/reads/ --output samplesheet.csv
```

Paired-end FASTQ files are detected and matched automatically using standard filename conventions (`_R1`/`_R2`, `_R1_001`/`_R2_001`, `_1`/`_2`). Files that do not match any paired pattern are treated as single-end.

### `cladebreaker prepare` options

| Option | Description |
|---|---|
| `--fastqs DIR` | Directory of FASTQ files (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) |
| `--assemblies DIR` | Directory of genome assemblies (`.fna`, `.fa`, `.fasta` and `.gz` variants) |
| `--output FILE` | Output samplesheet path (required) |
| `--recursive` | Search subdirectories recursively |
| `--prefix STR` | Only include files whose names begin with this prefix |

### Examples

```sh
# From paired-end reads
cladebreaker prepare --fastqs /data/reads/ --output samplesheet.csv

# From pre-assembled genomes
cladebreaker prepare --assemblies /data/assemblies/ --output samplesheet.csv

# Mixed reads and assemblies
cladebreaker prepare \
  --fastqs /data/reads/ \
  --assemblies /data/assemblies/ \
  --output samplesheet.csv

# Reads spread across subdirectories, filtered by prefix
cladebreaker prepare \
  --fastqs /data/reads/ \
  --recursive \
  --prefix SAMPLE_ \
  --output samplesheet.csv
```

---

## Step 3 — Run the Pipeline (`cladebreaker`)

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/database.pickle \
  --o \
  --topgenomes_count 5 \
  -profile conda
```

### Input/output options

| Option | Default | Description |
|---|---|---|
| `--input PATH` | required | Path to the input samplesheet CSV |
| `--outdir PATH` | `./results` | Directory where results will be saved |
| `--db PATH` | required | Path to the WhatsGNU database file |
| `--o` | — | Specify that the database is in ortholog mode |
| `--b` | — | Specify that the database is in basic mode |
| `--topgenomes_count INT` | `3` | Number of closest reference genomes to select per query |
| `--ref PATH` | — | Reference genome FASTA for reference-based alignment (enables Snippy mode) |
| `--email STR` | — | Email address for pipeline completion notification |

### Pangenome tool options

| Option | Description |
|---|---|
| `--pangenome_tool panaroo` | Use Panaroo for pangenome alignment (default) |
| `--pangenome_tool pirate` | Use PIRATE for pangenome alignment |
| `--pangenome_tool roary` | Use Roary for pangenome alignment (container profile required) |
| `--run_raxml` | Run RAxML-NG to build a maximum-likelihood phylogenetic tree |

### Assembly options (Shovill)

| Option | Default | Description |
|---|---|---|
| `--shovill_depth INT` | `150` | Target read depth for Shovill subsampling |
| `--shovill_gsize INT` | autodetect | Estimated genome size passed to Shovill |

### Annotation options (Prokka)

| Option | Description |
|---|---|
| `--proteins PATH` | Genbank file of trusted protein sequences for Prokka annotation |
| `--prodigal_tf PATH` | Pretrained Prodigal model file for gene prediction |

### Resource limit options

| Option | Default | Description |
|---|---|---|
| `--max_cpus INT` | `16` | Maximum CPUs for any single process |
| `--max_memory STR` | `128.GB` | Maximum memory for any single process |
| `--max_time STR` | `240.h` | Maximum wall time for any single process |

### Generic options

| Option | Default | Description |
|---|---|---|
| `--cleanup_workdir` | — | Delete the Nextflow `work/` directory after a successful run |
| `--force` | `false` | Overwrite existing output files |
| `--monochrome_logs` | — | Disable coloured log output |
| `--multiqc_config PATH` | — | Custom MultiQC config file |
| `--email_on_fail STR` | — | Email address for failure notifications only |

---

## Alignment Modes

### Pangenome alignment (default)

When no `--ref` is provided, Cladebreaker builds a core-genome alignment across all query samples and the downloaded reference genomes using your chosen pangenome tool. This approach:

- Detects the accessory genome as well as the core genome
- Does not require prior knowledge of a reference strain
- Works well for population-level analyses and clade assignment

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/db.pickle \
  --o \
  --topgenomes_count 5 \
  --pangenome_tool panaroo \
  --run_raxml \
  -profile conda
```

**Choosing a pangenome tool:**

| Tool | Notes |
|---|---|
| `panaroo` | Default. Robust graph-based approach with error correction. Works via conda. |
| `pirate` | Fast and scalable. Good for diverse or larger datasets. Works via conda. |
| `roary` | Classic rapid pangenome tool. **Requires a container profile** (Docker, Singularity, etc.) due to a Perl dependency conflict with conda. |

### Reference-based alignment (Snippy)

Add `--ref` to perform per-sample SNP calling against a reference genome with Snippy, followed by a core SNP alignment with Snippy-core. This approach:

- Is appropriate when a well-characterised reference genome is available
- Produces a compact SNP-only alignment ideal for closely related isolates
- Is commonly used for outbreak investigations

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/db.pickle \
  --o \
  --ref /path/to/reference.fna \
  --topgenomes_count 3 \
  --run_raxml \
  -profile conda
```

---

## Input Samplesheet Format

The input samplesheet is a four-column CSV file. Use `cladebreaker prepare` to generate it automatically, or create it manually:

```
sample,type,file_1,file_2
SAMPLE_PE,reads,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
SAMPLE_SE,reads,/path/to/reads.fastq.gz,
SAMPLE_ASM,assembly,/path/to/genome.fna,
```

| Column | Description |
|---|---|
| `sample` | Sample name. Used as the prefix for all output files. |
| `type` | Either `reads` (FASTQ input) or `assembly` (FASTA input). |
| `file_1` | Path to R1 FASTQ (paired/single-end) or FASTA assembly. |
| `file_2` | Path to R2 FASTQ for paired-end reads. Leave empty for single-end or assemblies. |

An example samplesheet is provided at `assets/samplesheet.csv`.

---

## Pipeline Overview

The following steps are run automatically for every execution:

```
Input validation
    └─ Samplesheet check → INPUT_CHECK

Quality control (reads only)
    └─ Raw read QC         → FastQC

Assembly (reads only)
    └─ De novo assembly    → Shovill

Assembly statistics
    └─ Assembly metrics    → assembly-scan

Genome annotation
    └─ Gene prediction     → Prokka

Closest-genome identification
    └─ Proteomic rarity    → WhatsGNU

Reference genome acquisition
    └─ Download from NCBI  → ncbi-genome-download + Prokka

Phylogenetic alignment (choose one mode)
    ├─ Pangenome mode      → Panaroo / PIRATE / Roary
    └─ Reference SNP mode  → Snippy + Snippy-core

Phylogenetic tree (optional)
    └─ ML inference        → RAxML-NG

Quality reporting
    └─ Aggregate report    → MultiQC
```

### Tool descriptions

| Tool | Role in pipeline |
|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Per-read quality metrics: base quality, adapter content, duplication |
| [Shovill](https://github.com/tseemann/shovill) | De novo assembly of Illumina short reads using SPAdes with read depth normalisation |
| [assembly-scan](https://github.com/rpetit3/assembly-scan) | Rapid summary statistics for genome assemblies (N50, contig count, total length) |
| [Prokka](https://github.com/tseemann/prokka) | Rapid prokaryotic genome annotation; generates GFF and FAA files consumed by downstream tools |
| [WhatsGNU](https://github.com/ahmedmagds/WhatsGNU) | Identifies the most closely related public genomes by scoring proteomic rarity across a compressed reference database |
| [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) | Downloads selected reference genome assemblies from NCBI FTP |
| [Panaroo](https://github.com/gtonkinhill/panaroo) | Graph-based pangenome pipeline; produces a core-gene multiple sequence alignment |
| [PIRATE](https://github.com/SionBayliss/PIRATE) | Fast, scalable pangenome tool; produces a core-gene alignment across diverged orthologues |
| [Roary](https://github.com/sanger-pathogens/Roary) | Classic rapid pangenome tool; produces a core-gene alignment (container only) |
| [Snippy](https://github.com/tseemann/snippy) | Reference-based SNP/indel calling from reads or assemblies; Snippy-core merges per-sample VCFs into a core SNP alignment |
| [RAxML-NG](https://github.com/amkozlov/raxml-ng) | Maximum-likelihood phylogenetic inference from core-gene or core-SNP alignments |
| [MultiQC](https://multiqc.info/) | Aggregates QC outputs from FastQC and other tools into a single HTML report |

---

## Output

All results are written under `--outdir`. The structure varies slightly depending on whether reads or assemblies were provided and which alignment mode was used.

### Per-sample outputs

| Directory | Contents |
|---|---|
| `{sample}/fastqc/` | FastQC HTML report and zip archive (reads input only) |
| `{sample}/shovill/` | De novo assembly FASTA and SPAdes log (reads input only) |
| `{sample}/assemblyscan/` | Assembly statistics JSON (contig count, N50, total length, GC%) |
| `{sample}/prokka/` | Genome annotation: GFF, FAA, GBK, and summary TSV |
| `{sample}/WhatsGNU/` | WhatsGNU rarity scores, top-genome report, and GCA accession list |

### Pangenome outputs (default mode)

| Directory | Contents |
|---|---|
| `panaroo/` | Pangenome results: core gene alignment (`core_gene_alignment.aln`), pan-genome graph, gene presence/absence matrix |
| `pirate/` | PIRATE pangenome results and core alignment (when `--pangenome_tool pirate`) |
| `roary/` | Roary pangenome results and core alignment (when `--pangenome_tool roary`) |

### Reference-based outputs (Snippy mode)

| Directory | Contents |
|---|---|
| `snippy/{sample}/` | Per-sample VCF, BED, and aligned FASTA |
| `snippy/core/` | Core SNP alignment (`core.full.aln`), VCF, and summary |

### Phylogenetic tree output (optional)

| Directory | Contents |
|---|---|
| `raxmlng/` | Best-scoring ML tree (`.raxml.bestTree`), bootstrap support trees, and log |

### Pipeline-wide outputs

| Directory | Contents |
|---|---|
| `multiqc/` | `multiqc_report.html`, per-tool parsed data, and static plot images |
| `pipeline_info/` | Nextflow execution report, timeline, trace file, and software versions YAML |

---

## Execution Profiles

Cladebreaker uses Nextflow profiles to configure the software environment. Specify one or more profiles with `-profile`:

| Profile | Description |
|---|---|
| `conda` | Install all dependencies via Conda/Mamba. Recommended for HPC without container support. |
| `docker` | Run each process in a Docker container. Fully reproducible. |
| `singularity` | Run each process in a Singularity container. Preferred on shared HPC systems. |
| `podman` | Run each process with Podman (rootless containers). |
| `shifter` | Run each process with Shifter (NERSC and similar systems). |
| `charliecloud` | Run each process with Charliecloud. |
| `local` | Laptop/workstation preset: caps resources at 8 CPUs, 14 GB RAM, 48 h. Use with `conda,local`. |
| `test` | Minimal test dataset for pipeline validation. Requires no additional parameters. |

Multiple profiles can be combined:

```sh
# Conda on a laptop
cladebreaker --input samplesheet.csv --db db.pickle --o -profile conda,local

# Singularity on an HPC cluster
cladebreaker --input samplesheet.csv --db db.pickle --o -profile singularity -c cluster.config

# Docker with a custom resource config
cladebreaker --input samplesheet.csv --db db.pickle --o -profile docker -c custom.config
```

> **Note on Roary**: Roary cannot be installed via conda due to a Perl dependency conflict. Use `docker`, `singularity`, or another container profile when `--pangenome_tool roary` is set.

---

## Configuration

### Running on an HPC cluster

Supply a Nextflow config file with your cluster settings via `-c`. A typical SLURM config looks like:

```nextflow
process {
    executor = 'slurm'
    queue    = 'short'
    clusterOptions = '--account=myproject'
}
```

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/db.pickle \
  --o \
  -profile conda \
  -c cluster.config
```

### Adjusting resource limits

Cap the maximum resources any single process can request using the pipeline's resource parameters:

```sh
cladebreaker --input samplesheet.csv --db db.pickle --o \
  --max_cpus 32 \
  --max_memory 256.GB \
  --max_time 120.h \
  -profile conda
```

### Overriding per-process resources

To tune resources for a specific process, add a `withName` block to a custom config:

```nextflow
process {
    withName: PANAROO {
        cpus   = 24
        memory = 64.GB
        time   = 24.h
    }
}
```

### Resuming a run

Nextflow caches all completed process outputs in the `work/` directory. If a run fails or is interrupted, resume from where it left off with `-resume`:

```sh
cladebreaker --input samplesheet.csv --db db.pickle --o -profile conda -resume
```

### Running in the background

Use the Nextflow `-bg` flag or a terminal multiplexer (`screen`, `tmux`) to detach from the session:

```sh
cladebreaker --input samplesheet.csv --db db.pickle --o -profile conda -bg
```

### Keeping memory use low

Nextflow's JVM can accumulate memory over long runs. Add the following to your `~/.bashrc` to cap it:

```sh
export NXF_OPTS='-Xms1g -Xmx4g'
```

---

## Utility Subcommands

### `cladebreaker citations`

Print citations for all tools used by Cladebreaker:

```sh
cladebreaker citations
```

### `cladebreaker version`

Print the installed Cladebreaker version:

```sh
cladebreaker version
# cladebreaker v0.3.1
```

---

## Citations

If you use Cladebreaker in your research, please cite:

> Feder AF, Planet PJ (2022). cladebreaker: a Nextflow pipeline for whole-genome
> sequencing-based clade assignment of bacterial isolates.
> https://github.com/andriesfeder/cladebreaker

For a full list of per-tool citations, run `cladebreaker citations`. The underlying tools are the foundation of this pipeline — please also cite them when publishing:

| Tool | Reference |
|---|---|
| Nextflow | Di Tommaso *et al.* (2017) *Nat Biotechnol* [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| nf-core | Ewels *et al.* (2020) *Nat Biotechnol* [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x) |
| WhatsGNU | Moustafa & Planet (2020) *Genome Biol* [doi:10.1186/s13059-020-01965-w](https://doi.org/10.1186/s13059-020-01965-w) |
| Prokka | Seemann (2014) *Bioinformatics* [doi:10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153) |
| Panaroo | Tonkin-Hill *et al.* (2020) *Genome Biol* [doi:10.1186/s13059-020-02090-4](https://doi.org/10.1186/s13059-020-02090-4) |
| PIRATE | Bayliss *et al.* (2019) *Gigascience* [doi:10.1093/gigascience/giz119](https://doi.org/10.1093/gigascience/giz119) |
| Roary | Page *et al.* (2015) *Bioinformatics* [doi:10.1093/bioinformatics/btv421](https://doi.org/10.1093/bioinformatics/btv421) |
| RAxML-NG | Kozlov *et al.* (2019) *Bioinformatics* [doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305) |
| Shovill | Seemann (GitHub) [tseemann/shovill](https://github.com/tseemann/shovill) |
| Snippy | Seemann (GitHub) [tseemann/snippy](https://github.com/tseemann/snippy) |
| FastQC | Andrews (Babraham) [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| MultiQC | Ewels *et al.* (2016) *Bioinformatics* [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |
| ncbi-genome-download | Blin *et al.* (GitHub) [kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) |
| assembly-scan | Petit (GitHub) [rpetit3/assembly-scan](https://github.com/rpetit3/assembly-scan) |

---

## Support

- **Bug reports and feature requests**: please open an issue at [github.com/andriesfeder/cladebreaker/issues](https://github.com/andriesfeder/cladebreaker/issues)
- **Questions**: check existing issues before opening a new one — your question may already have an answer

---

## License

Cladebreaker is released under the [MIT License](LICENSE).
