# cladebreaker: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`

- **`bin/monophyly.py`**: Reports whether each group is a clade, the fraction of an ensemble in which it is, the number of separate clusters a broken group falls into, and the named clade-breaking tips
- **`bin/rosenberg.py`**: Rosenberg (2007) `P_A` and `P_AB` for chance monophyly, Holm-adjusted across groups, plus the k-group joint probability of Zhu, Degnan & Steel (2011) Theorem 5.1 evaluated by an O(3^k) subset recursion rather than enumerating (2k-3)!! trees. Computed in log space so realistic sample sizes neither overflow the binomial nor underflow the float
- **`bin/slatkin_maddison.py`**: Slatkin & Maddison (1989) migration test by Fitch parsimony, with the polytomy and missing-data generalisations, permutation p-value, and the null mean/SD for context. Independent of rooting
- **`bin/snp_separation.py`**: Within- and between-group SNP distance summaries with `ratio_of_means` and `gap` (min-between over max-within), a permutation p-value, and a pairwise between-group table. Ratios are defined by this pipeline, not drawn from the literature, and say so in every output
- **`bin/analyze_pdf.py`** and the `ANALYZE_PDF` process: Renders the results as a PDF -- the rooted tree the statistics were computed on with tips marked by group, then one section per test, then the SNP distance ranges and caveats (`conda-forge::matplotlib-base`). Colour is used for at most three groups because a fourth categorical slot falls below the normal-vision separation floor; beyond that identity moves to marker shape, and tip labels always carry it regardless
- **`--no-report` parameter**: Skips the PDF (`--no_report` is accepted too). The PDF is also skipped automatically when `--tests` does not select `monophyly`, which is the test that writes the rooted tree the figure draws
- **`monophyly.py`**: Also writes `<prefix>.rooted.nwk`, the tree the tests ran on, so the report figure draws that rather than the input file
- **`bin/analyze_report.py`**: Applies the decision path — Rosenberg where the groups are clades, Slatkin-Maddison where they are not — and writes a combined verdict. All tests run regardless; the decision only chooses which result leads
- **`bin/cladebreaker_phylo.py`**: Shared tree loading, rooting, groups parsing, tree indexing, monophyly and label-permutation driver, extracted from `gsi.py` so every analyze tool builds on the same primitives
- **`--tests` parameter**: `auto` (the default, runs the whole decision path) or a comma-separated subset of `gsi,monophyly,rosenberg,slatkin_maddison,snp_separation`
- **`--run_analysis` parameter**: Runs the full decision path on the tree built by the main pipeline, superseding `--run_gsi` which runs the gsi alone
- **`--alignment` / `--distances` parameters**: Core alignment (via `SNP_DISTS`, `bioconda::snp-dists=1.2.0`) or a precomputed matrix, for the SNP separation test
- **`--gsi_alpha` and `--gsi_max_joint_groups` parameters**: Significance level for the report, and the group-count cap on the joint monophyly probability
- **`modules/local/analyze/`**: `MONOPHYLY`, `ROSENBERG`, `SLATKIN_MADDISON`, `SNP_DISTS`, `SNP_SEPARATION` and `ANALYZE_REPORT` processes, publishing to `<outdir>/analyze/`
- **`ANALYZE` workflow** (`workflows/analyze.nf`): New entry point that runs statistics on an existing phylogenetic tree without rebuilding it — `nextflow run main.nf --workflow ANALYZE --tree <tree.nwk> --groups <groups.tsv>`
- **`bin/gsi.py`**: Genealogical sorting index (gsi) of [Cummings, Neel & Shaw (2008)](https://doi.org/10.1111/j.1558-5646.2008.00442.x), quantifying how far each labeled group on a rooted tree has sorted towards exclusive ancestry (0 = mixed, 1 = monophyletic), with a label-permutation p-value and the ensemble statistic `gsi_T` across bootstrap or posterior trees. Outputs TSV, JSON, and a MultiQC custom-content table
- **`cladebreaker analyze` subcommand** (`bin/cladebreaker/cladebreaker`): Dispatches to the `ANALYZE` workflow, alongside the existing `build` and `prepare` subcommands, with its own `--help`
- **`bin/gsi_groups.py`**: Builds the tip-to-group table by resolving tree tip labels against the pipeline's query and reference sample IDs, tolerating the suffixes annotation and pangenome steps append while refusing ambiguous matches
- **`modules/local/gsi/main.nf`** and **`modules/local/gsi/groups.nf`**: New `GSI` and `GSI_GROUPS` processes (`bioconda::dendropy=5.0.13`, container `quay.io/biocontainers/dendropy:5.0.13--pyhdfd78af_0`)
- **`--run_gsi` parameter**: Runs the gsi on the tree built by the main pipeline, scoring the input isolates against the reference genomes WhatsGNU selected. Requires `--run_raxml`
- **gsi parameters**: `--tree`, `--groups`, `--gsi_root` (`as-is`/`midpoint`/`outgroup`), `--gsi_outgroup`, `--gsi_permutations`, `--gsi_seed`, `--gsi_groups_to_test`, `--gsi_ignore_unlabeled`, `--gsi_query_label`, `--gsi_reference_label`
- **`--workflow` parameter**: Now declared in `nextflow.config` and the schema (`BUILD` or `ANALYZE`); it was previously only ever set from the command line
- **`modules/nf-core/modules/raxmlng/main.nf`**: Added optional `bootstrap_trees` output (`*.raxml.bootstraps`) so bootstrap replicates can feed the ensemble gsi
- **`conf/test_analyze.config`** and the `test_analyze` profile: Runs the whole decision path against the worked example from Fig. 1 of Cummings, Neel & Shaw (2008). One group there is a clade and the other is not, so a single run exercises both branches -- Rosenberg for the first, Slatkin-Maddison for the second
- **`tests/gsi/cummings_fig1.aln`**: A 30-column alignment over the same tips, built so every distance is checkable by hand (2 SNPs within group a, 4 within group b, 9 between), giving the test profile coverage of the SNP separation test
- **`tests/gsi/`**: Test suite for all the analyze scripts (75 tests), including Rosenberg's 28 published table values and four worked examples, the Zhu et al. `p(2,2,2)=2/225` example with its reduction to `P_AB` at k=2, Fitch parsimony cross-checked against an exhaustive minimum over every internal-node labeling, rooting-invariance of the Slatkin-Maddison statistic checked at every edge, the gsi paper's Fig. 1 and Fig. 2 worked examples with an exhaustive enumeration of the null distribution, midpoint-rooting correctness, and the sample-size power floors that limit all three permutation tests
- **`modules/local/bakta/main.nf`**: New `BAKTA` process for genome annotation using [Bakta](https://github.com/oschwengers/bakta) as an alternative to Prokka in both the main pipeline and `cladebreaker build` (`bioconda::bakta`, container `quay.io/biocontainers/bakta:1.9.4`)
- **`--annotator` parameter**: Selects the annotation tool — `prokka` (default) or `bakta` — for both user assemblies and downloaded reference genomes in the main pipeline, and for reference genomes in the `cladebreaker build` workflow
- **`--bakta_db PATH` parameter**: Path to a local Bakta database directory; required when `--annotator bakta` is set. Database can be downloaded with `bakta_db download`

### `Changed`

- **`bin/gsi.py`**: Implements midpoint rooting directly rather than calling DendroPy's `Tree.reroot_at_midpoint()`. As of DendroPy 5.0.13 that method walks up from only one end of the diameter and raises a bare `AssertionError` when floating-point drift leaves distance remaining at the MRCA; across 400 random trees it failed on 94 and returned an off-midpoint rooting on a further 46
- **`bin/gsi.py`**: Now imports its tree handling from the new `bin/cladebreaker_phylo.py` rather than defining it. Verified behaviour-preserving: the existing tests pass and output on a real 77-tip tree is byte-identical, seeded permutation p-values included
- **`nextflow.config`**: Added `run_gsi`, `gsi_ignore_unlabeled` and `run_analysis` to `schema_ignore_params`, alongside the existing `run_raxml`. Boolean flags passed bare on the command line arrive as the string `"true"`, which the schema validator rejects as a `Boolean`; ignoring them skips validation only, and they remain documented in `--help`

### `Fixed`

- **`main.nf`**: Replaced `-entry CLADEBREAKER_BUILD` workflow selection with a `params.workflow` branch inside the entry `workflow {}` block — the `-entry` option is not supported by the Nextflow 26.x strict parser
- **`bin/cladebreaker-build.py`**: Updated `launch_build_workflow()` to pass `--workflow BUILD` instead of `-entry CLADEBREAKER_BUILD` to match the new entry-point pattern
- **`modules/local/build/fetch_accessions.nf`**: Replaced bash `cat <<-END_VERSIONS` heredoc with Python `open()` to write `versions.yml` — the heredoc caused a `SyntaxError` because the entire script block is executed as Python (shebang `#!/usr/bin/env python3`)
- **`modules/local/build/fetch_accessions.nf`**: Updated conda directive from `bioconda::biopython=1.79` to `conda-forge::biopython` — biopython is distributed via `conda-forge`, not `bioconda`, and version `1.79` was not resolvable

## v0.3.0 - 2026-05-29

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
