# nf-core/cladebreaker: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/cladebreaker, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Added `local` profile (`max_memory = 14.GB`, `max_cpus = 8`, `max_time = 48.h`) for running the pipeline on a laptop/workstation with `-profile conda,local` ([#](https://github.com/nf-core/cladebreaker))

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
- **`modules/nf-core/modules/prokka/main.nf`**: Updated conda directive from `prokka=1.14.5` to `prokka=1.14.6` — the 1.14.5 package has an unresolvable Perl dependency conflict in current bioconda

### `Dependencies`

- Updated Prokka: `1.14.5` → `1.14.6`

### `Deprecated`
