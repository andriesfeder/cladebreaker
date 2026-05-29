# Cladebreaker: Phylogenetic Analysis Pipeline

Cladebreaker is a comprehensive pipeline for the phylogenetic analysis of bacterial genomes. It identifies closely related genomes using [WhatsGNU](https://github.com/ahmedmagds/WhatsGNU), creates pangenome alignments with [Panaroo](https://github.com/gtonkinhill/panaroo) (default), [PIRATE](https://github.com/SionBayliss/PIRATE), or [Roary](https://github.com/sanger-pathogens/Roary) (container only), and supports reference-based alignments via Snippy — enabling genomic comparisons and evolutionary studies for any bacterial species.

## Installation

Cladebreaker can be installed using Mamba (recommended) or Conda:

```sh
mamba create -n cladebreaker -c bioconda cladebreaker
mamba activate cladebreaker
```

Verify the installation:

```sh
cladebreaker version
```

## Usage

Cladebreaker is invoked through a single entry point with several subcommands:

```
cladebreaker [options]            Run the main analysis pipeline
cladebreaker build [options]      Download or build a WhatsGNU database
cladebreaker prepare [options]    Generate an input samplesheet
cladebreaker citations            Print citations for all tools
cladebreaker version              Print the cladebreaker version
```

---

## Step 1 — Prepare a database (`cladebreaker build`)

Cladebreaker requires a WhatsGNU database to identify closely related genomes.
The `build` subcommand handles this automatically: it checks whether a pre-built
database is available for your species, and if so prompts you to download it.
If no pre-built database exists, it builds one from NCBI reference genomes.

```sh
cladebreaker build --taxid 1280
```

Pre-built databases are available for the following species:

| Species | TaxID | Modes available |
|---|---|---|
| *Staphylococcus aureus* | 1280 | ortholog |
| *Staphylococcus epidermidis* | 1282 | basic |
| *Klebsiella pneumoniae* | 573 | ortholog, basic |
| *Pseudomonas aeruginosa* | 287 | ortholog, basic |
| *Escherichia coli* | 562 | basic |
| *Salmonella enterica* | 28901 | basic |
| *Mycobacterium tuberculosis* | 1773 | ortholog |
| *Clostridioides difficile* | 1496 | ortholog |

For any species not listed above, `cladebreaker build` will automatically
download reference genomes from NCBI and construct a database.

**Key options:**

| Option | Description |
|---|---|
| `--taxid INT` | NCBI TaxID of the target species (required) |
| `--ortholog` | Build/download an ortholog-mode database (default: basic) |
| `--genome_count INT` | Max genomes for a custom build (default: all available) |
| `--outdir PATH` | Where to save the database (default: `./cladebreaker_databases`) |
| `-profile STR` | Nextflow profile for a custom build (e.g. `conda`, `docker`) |
| `-c PATH` | Additional Nextflow config file |

```sh
# Download a pre-built basic database for K. pneumoniae
cladebreaker build --taxid 573

# Build a custom database for A. baumannii using 200 genomes
cladebreaker build --taxid 470 --genome_count 200 --outdir /data/databases -profile conda

# Build an ortholog database for S. aureus
cladebreaker build --taxid 1280 --ortholog
```

The `build` subcommand will display a confirmation summary before starting any
download or build, including the species name, genome count, database mode, and
output location.

---

## Step 2 — Prepare a samplesheet (`cladebreaker prepare`)

Cladebreaker requires a CSV samplesheet describing your input samples.
The `prepare` subcommand generates one automatically from a directory of
FASTQ files and/or genome assemblies.

```sh
cladebreaker prepare --fastqs /path/to/reads/ --output samplesheet.csv
```

Paired-end FASTQ files are detected automatically using standard R1/R2 filename
conventions (`_R1`/`_R2`, `_R1_001`/`_R2_001`, `_1`/`_2`).

**Key options:**

| Option | Description |
|---|---|
| `--fastqs DIR` | Directory of FASTQ files (`.fastq.gz`, `.fq.gz`) |
| `--assemblies DIR` | Directory of genome assemblies (`.fna`, `.fa`, `.fasta`) |
| `--output FILE` | Output samplesheet path (required) |
| `--recursive` | Search subdirectories recursively |
| `--prefix STR` | Only include files whose names begin with this prefix |

```sh
# From a directory of paired-end reads
cladebreaker prepare --fastqs /data/reads/ --output samplesheet.csv

# From a directory of assemblies
cladebreaker prepare --assemblies /data/assemblies/ --output samplesheet.csv

# Mixed reads and assemblies
cladebreaker prepare --fastqs /data/reads/ --assemblies /data/assemblies/ --output samplesheet.csv
```

The generated samplesheet can be passed directly to `cladebreaker` with `--input`.

---

## Step 3 — Run the pipeline

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir /PATH/TO/results \
  --db /PATH/TO/database.pickle \
  --o \
  --topgenomes_count 5 \
  -profile conda \
  -c /PATH/TO/cluster.config
```

**Key options:**

| Option | Description |
|---|---|
| `--input PATH` | Path to the input samplesheet CSV (required) |
| `--outdir PATH` | Directory where results will be saved |
| `--db PATH` | Path to the WhatsGNU database pickle file (required) |
| `--o` | Use ortholog database mode |
| `--b` | Use basic database mode |
| `--topgenomes_count INT` | Number of top genomes to select per query (default: 3) |
| `--ref PATH` | Path to a reference genome for reference-based alignment |
| `--run_raxml` | Run RAxML-NG to generate a phylogenetic tree |
| `--pangenome_tool STR` | Pan-genome tool: `panaroo` (default), `pirate`, or `roary` (container only) |
| `--coverage INT` | Coverage threshold for genome filtering |
| `--force` | Overwrite existing output files |
| `-profile STR` | Nextflow execution profile (e.g. `conda`, `docker`, `singularity`) |
| `-c PATH` | Additional Nextflow config file (e.g. cluster configuration) |
| `--cleanup_workdir` | Remove the Nextflow work directory after a successful run |

### Pangenome alignment (default)

By default, Cladebreaker builds a pangenome alignment using Panaroo. Use `--pangenome_tool pirate` to switch to PIRATE, or `--pangenome_tool roary` to use Roary (requires a container profile — Roary cannot be installed via conda due to a Perl dependency conflict):

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/db.pickle \
  --o \
  --topgenomes_count 5 \
  -profile conda
```

### Reference-based alignment

Add `--ref` to perform a reference-based SNP alignment with Snippy instead:

```sh
cladebreaker \
  --input samplesheet.csv \
  --outdir results/ \
  --db /path/to/db.pickle \
  --o \
  --ref /path/to/reference.fna \
  -profile conda
```

---

## Input file format

The input samplesheet is a CSV file with the following columns:

```
sample,type,file_1,file_2
SAMPLE_PE,reads,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
SAMPLE_SE,reads,/path/to/R1.fastq.gz,
SAMPLE_ASM,assembly,/path/to/genome.fna,
```

Use `cladebreaker prepare` to generate this file automatically.
An example samplesheet is available at `assets/samplesheet.csv`.

---

## Output

Cladebreaker produces the following results under `--outdir`:

| Output | Description |
|---|---|
| `{sample}/WhatsGNU/` | WhatsGNU reports and top-genome lists |
| `{sample}/prokka/` | Prokka genome annotations |
| `{sample}/shovill/` | Shovill assemblies (if FASTQ input) |
| `{sample}/fastqc/` | FastQC read quality reports |
| `panaroo/` | Panaroo pangenome results and core gene alignment (default) |
| `pirate/` | PIRATE pangenome results (when `--pangenome_tool pirate`) |
| `roary/` | Roary pangenome results (when `--pangenome_tool roary`) |
| `snippy/` | Snippy per-sample variant calls and core SNP alignment |
| `raxmlng/` | RAxML-NG phylogenetic tree (if `--run_raxml`) |
| `multiqc/` | MultiQC summary report |
| `pipeline_info/` | Nextflow execution timeline, trace, and reports |

---

## Citations

For a full list of citations, run:

```sh
cladebreaker citations
```

If you use Cladebreaker in your research, please cite:

> Feder AF, Planet PJ (2022). cladebreaker: a Nextflow pipeline for whole-genome
> sequencing-based clade assignment of bacterial isolates.
> https://github.com/andriesfeder/cladebreaker

Please also cite the underlying tools: WhatsGNU, Nextflow, Prokka, Panaroo, PIRATE,
Roary, Snippy, RAxML-NG, Shovill, FastQC, MultiQC, and ncbi-genome-download.
Run `cladebreaker citations` for the full reference list.

---

## Support

Please open an issue at https://github.com/andriesfeder/cladebreaker/issues
for bug reports or feature requests.

## License

Cladebreaker is released under the MIT License. See the `LICENSE` file for details.
