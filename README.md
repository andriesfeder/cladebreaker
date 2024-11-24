# Cladebreaker: Phylogenetic Analysis Pipeline

Cladebreaker is a comprehensive, community-driven pipeline designed for the phylogenetic analysis of bacterial genomes. This tool helps identify closely related genomes, create pangenome alignments, and perform reference-based alignments to support genomic comparisons and evolutionary studies.

## Installation

Cladebreaker can be installed using either Conda or Mamba. We recommend using Mamba for faster dependency resolution, but Conda will work as well.

```sh
mamba create -n cladebreaker -c bioconda cladebreaker
mamba activate cladebreaker
```

To verify the installation, check the version (the latest version is 0.2.3):

```sh
cladebreaker --version
```

## Requirements

To run Cladebreaker successfully, you will need a WhatsGNU database, which can be downloaded from the [WhatsGNU GitHub page](https://github.com/ahmedmagds/WhatsGNU). For Staphylococcus aureus analysis, we recommend using the Staph. aureus ortholog database as it enables additional downstream analysis. Alternatively, you can use the basic WhatsGNU database, but it will require slight modifications to the commands below.

## Running Cladebreaker

The following command runs Cladebreaker to identify the top 5 most closely related genomes to each query genome and create a pangenome alignment using Roary:

```sh
cladebreaker \
  --input /PATH/TO/cladebreaker_input.csv \
  --outdir /PATH/TO/cladebreaker_results \
  --coverage 100 \
  --db /PATH/TO/WhatsGNU_Sau_Ortholog/Sau_Ortholog_10350.pickle \
  --o \
  --topgenomes_count 5 \
  --force \
  -profile conda \
  -c /PATH/TO/penn_cluster.config
```

- **--input**: Path to the input CSV file containing metadata for each query genome.
- **--outdir**: Directory where the results will be saved.
- **--coverage**: Coverage threshold for genome selection.
- **--db**: Path to the WhatsGNU database.
- **--o**: Use ortholog-based database (replace with `--b` for basic database).
- **--topgenomes_count**: Number of top genomes to select for each query.
- **--force**: Force overwrite of the output directory if it already exists.
- **--profile conda**: Use the Conda environment profile for reproducibility.
- **-c**: Path to the cluster configuration file.

### Reference-Based Alignment

If you prefer to perform a reference-based alignment instead of a pangenome analysis, simply add the `--ref` flag with the path to the reference genome you'd like to use:

```sh
--ref /Documents/cladebreaker/N315.fna
```

## Input File Format

The input for Cladebreaker should be in CSV format, containing information about each query genome. An example input CSV file is included in this repository to help you get started.

## Output

Cladebreaker produces the following results:
- **Top Related Genomes**: Identification of the most closely related genomes to each query.
- **Pangenome Alignment**: Creation of a pangenome alignment using Roary.
- **Optional Reference-Based Alignment**: Alignment against a provided reference genome if specified.

## Citation

If you use Cladebreaker in your research, please cite the following publications:
- Ewels PA, et al. "The nf-core framework for community-curated bioinformatics pipelines." *Nat Biotechnol*. 2020.
- Di Tommaso P, et al. "Nextflow enables reproducible computational workflows." *Nat Biotechnol*. 2017.

For more detailed citation information, refer to the [CITATIONS.md](https://github.com/andriesfeder/cladebreaker/blob/master/CITATIONS.md) file.

## Support and Contributions

We welcome contributions from the community! Please fork the repository and submit a pull request, or open an issue for any suggestions or bug reports. For further assistance, refer to the [nf-core Slack](https://nfcore.slack.com/channels/cladebreaker) channel.

## License

Cladebreaker is released under the MIT License. See the `LICENSE` file for more information.

