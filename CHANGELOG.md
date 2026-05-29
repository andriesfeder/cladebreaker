# cladebreaker: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of cladebreaker, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `local` profile (`max_memory = 14.GB`, `max_cpus = 8`, `max_time = 48.h`) for running the pipeline on a laptop/workstation with `-profile conda,local` ([#](https://github.com/andriesfeder/cladebreaker))

### `Fixed`

#### Nextflow 26.x Compatibility

- **`nextflow.config`**: Removed `check_max()` function definition (method-style `def` banned in NF 26.x config); replaced `try/catch` around `includeConfig` with a ternary conditional (try-catch banned in config); inlined `trace_timestamp` variable directly into file paths (top-level variable declarations banned); skipped remote HTTP `includeConfig` to avoid SSL errors when running offline; added `conda.enabled = true` and `conda.useMamba = true` to the conda profile
- **`conf/base.config`**: Replaced all `check_max(obj, type)` calls with inline `[val * task.attempt, params.max_x as Type].min()` expressions — lib Groovy classes are inaccessible from config file closures at runtime in NF 26.x
- **`main.nf`**: Removed top-level statements (banned in NF 26.x); moved `WorkflowMain.initialise()` into the `workflow {}` entry block
- **`workflows/cladebreaker.nf`**: Moved all validation and initialization code inside `workflow CLADEBREAKER {}`; explicitly captured `workflow`, `params`, `projectDir`, and `log` as `def` locals for use in `onComplete` closure (implicit variables not reliably available inside named-workflow closures in NF 26.x)
- **`subworkflows/local/input_check.nf`** and **`gather_genomes.nf`**: Removed all `import` statements (banned in NF 26.x `.nf` files); replaced with fully-qualified class names; changed all typed local variable declarations (`String`, `Path`, `File`, `final File`) to `def`; replaced `for` loops with `.each {}` closures (for loops banned in NF 26.x)
- **`modules/local/cladebreaker/qc_reads.nf`**: Removed top-level `RESOURCES` and `options` variable declarations; wrapped `publishDir` path in closure; changed `path()` to `file()` in script context
- **`modules/local/snippy/snippycore.nf`**: Replaced `for(i in paths)` loop with `paths instanceof List ? paths.join(' ') : paths`
- **`modules/nf-core/modules/ncbigenomedownload/main.nf`**: Moved `errorStrategy` directive before the `input:` block (was after `output:`, where it was parsed as an output qualifier); fixed syntax from `in 1` to `in [1]`
- **All modules with `${meta.id}` in `publishDir`** (assemblyscan, bbmap/bbduk, fastqc, prokka, shovill, snippy, whatsgnu): Wrapped `publishDir` paths in closures for lazy evaluation
- **`lib/NfcoreTemplate.groovy`**: Fixed unescaped backslashes in ASCII art logo (caused Groovy parse errors)
- **`lib/nf/functions.nf`**: Added `def` to undeclared `ext` variable; removed unused `get_resources()` function which referenced undefined internal methods
- **`nextflow_schema.json`**: Changed `o`, `b`, and `force` parameter types to `["boolean", "string"]` and `topgenomes_count` to `["integer", "string"]` to support CLI flag passing; extended `schema_ignore_params` to cover all pipeline-specific params

#### Runtime / Environment Fixes

- **`conda/meta.yaml`**: Replaced `conda` run dependency with `biopython` — having `conda` as a package dependency installed an isolated conda binary inside the environment that returned null for `conda_prefix`, causing Nextflow's conda activation scripts to fail with `/bin/activate: No such file or directory`
- **`nextflow.config`**: Added `PERL5LIB` to the `env {}` block to fix a macOS conda issue where the Perl binary has a compiled-in `@INC` path that does not include where conda actually installs modules; the value uses `${CONDA_PREFIX}` which bash expands at task runtime after conda activation
- **`nextflow.config`**: Added `PATH` to the `env {}` block to fix `ModuleNotFoundError: No module named 'yaml'` in CUSTOM_DUMPSOFTWAREVERSIONS (and analogous tool-shadowing bugs) when Nextflow is launched from within an active conda environment. Root cause: Nextflow's `source activate` in each task does a *stacked* activation that appends the task env's `bin/` instead of prepending it, so executables from the outer env (e.g. the user's `cladebreaker` env) shadow the task env's tools. `CONDA_PREFIX` is set correctly to the task env by conda even after a stacked activation; re-prepending `${CONDA_PREFIX}/bin` in the `env {}` block (which runs after activation) restores the correct lookup order for all tasks
- **`modules/nf-core/modules/prokka/main.nf`**: Updated conda directive from `prokka=1.14.5` to `prokka=1.14.6` — the 1.14.5 package has an unresolvable Perl dependency conflict in current bioconda
- **`modules/nf-core/modules/multiqc/main.nf`**: Upgraded from `multiqc=1.12` to `multiqc=1.35` — versions prior to 1.21 import `pkg_resources` (from `setuptools`), which is no longer included by default when conda resolves Python 3.12+ environments; multiqc 1.21+ replaced `pkg_resources` with `importlib.metadata`
- **`modules/nf-core/modules/panaroo/main.nf`**: Pinned `conda-forge::biopython<=1.86` alongside `panaroo=1.6.0` — BioPython 1.87 changed `FastaIterator` to read the first line strictly and reject it if it is not `>`, breaking panaroo's GFF FASTA extraction which splits on `##FASTA` and leaves a leading newline before the first record; BioPython ≤1.86 handled leading whitespace correctly
- **`modules/nf-core/modules/raxmlng/main.nf`**: Upgraded from `raxml-ng=1.0.3` to `raxml-ng=2.0.1`; updated container tags to `2.0.1--haec14ce_0`
- **`main.nf`**: Removed redundant `NFCORE_CLADEBREAKER` and `NFCORE_BUILD` wrapper workflows; entry-point `workflow {}` block now calls `CLADEBREAKER()` directly; fixed `BUILD` include to match the actual workflow name in `build.nf`
- **`workflows/build.nf`**: Renamed `workflow CLADEBREAKER_BUILD` to `workflow BUILD` for a cleaner entry-point name (`-entry BUILD`)
- **`assets/email_template.html`**, **`nextflow.config`**: Rebranded pipeline references from `nf-core/cladebreaker` to `cladebreaker`
- **Pan-genome tool selection** (`modules/nf-core/modules/panaroo/main.nf`, `modules/nf-core/modules/pirate/main.nf`, `modules/nf-core/modules/roary/main.nf`, `workflows/cladebreaker.nf`, `nextflow.config`): Replaced hard-coded Roary with a configurable pan-genome tool selection via `--pangenome_tool`. Added Panaroo 1.6.0 as the new default (resolves the bioconda Perl dependency conflict that prevents Roary from installing under conda). PIRATE 1.0.5 is available as an alternative (`--pangenome_tool pirate`). Roary is retained for container-based runs (`--pangenome_tool roary`); a startup warning is emitted when Roary is selected with `-profile conda`. All three tools now use the same GFF input format (collected files staged in the work directory). Updated PIRATE to 1.0.5 and added `publishDir` to PIRATE and Roary modules
- **`lib/WorkflowCladebreaker.groovy`**: Added startup warning when `--pangenome_tool roary` is combined with `-profile conda`

### `Dependencies`

- Updated Prokka: `1.14.5` → `1.14.6`
- Updated MultiQC: `1.12` → `1.35`
- Updated RAxML-NG: `1.0.3` → `2.0.1`

### `Deprecated`
