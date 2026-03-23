# miND - miRNA NGS Discovery Pipeline

miND is a small RNA-seq analysis pipeline for microRNA biomarker discovery studies. It processes raw sequencing data through quality control, hierarchical read mapping, miRNA quantification (via miRDeep2), classification of other small RNA species, and differential expression analysis, producing a comprehensive interactive HTML report.

## Features

- **End-to-end workflow**: From raw FASTQ/BAM files to differential expression results in a single pipeline run
- **Interactive HTML report**: Structured report with QC metrics, RNA class distributions, heatmaps, PCA/t-SNE, and DE results, including export for publication-ready figures
- **Excel-based contrast sheet**: Define sample groups and statistical comparisons in a familiar spreadsheet format
- **miRNA quantification**: Based on miRDeep2 with isomiR-level resolution
- **Small RNA classification**: Quantification of tRNAs, piRNAs, rRNAs, snRNAs, snoRNAs, lncRNAs, and other ncRNAs via RNAcentral
- **Differential expression**: edgeR with an adapted independent filtering method
- **Reproducible**: Snakemake workflow with conda-managed dependencies

## Requirements

- Linux (developed and tested on Debian)
- [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
- Hardware: 4+ CPU cores, 8 GB RAM (scales with available resources)

## Installation

A step-by-step installation and setup guide, including the preparation of the reference data repository, is available as a protocol on protocols.io:

**[miND pipeline AWS EC2 installation and setup V.2](https://dx.doi.org/10.17504/protocols.io.b3f6qjre)**

The protocol describes setup on an AWS EC2 instance but applies to any Linux system. Only OS-specific commands (e.g., package installation via `apt`) need to be adapted for other distributions.

### Quick overview

```bash
# Clone the repository
git clone https://github.com/TAmiRNA/miND.git
cd miND

# Install conda environments (handled automatically by Snakemake)
# Build the reference data repository (run once)
bash repository/build.sh

# The pipeline is now ready to use
```

## Usage

### 1. Prepare the contrast sheet

Copy `SampleContrastSheet.example.xlsx` and fill in three sheets:

1. **Project details**: Project name, species, adapter sequence, significance cutoffs
2. **Sample group matrix**: Sample names and group assignments (up to 5 grouping variables)
3. **Contrast selection**: Select group comparisons for differential expression analysis

### 2. Configure and run

Edit `config.yaml` to set paths to your data and contrast sheet, then run:

```bash
bash run.sh
```

The pipeline will process all samples and generate the HTML report in the output directory.

## Output

- **HTML report**: Interactive report with all QC, exploration, and DE results
- **Count matrices**: Raw and RPM-normalized miRNA read counts (CSV)
- **Mapping statistics**: Per-sample RNA class distributions
- **DE results**: Tables of differentially expressed miRNAs for each contrast
- **miRDeep2 output**: Novel miRNA predictions and isomiR mappings

An example report generated from a public dataset is available in the [examples](examples/) directory.

## Citation

If you use miND in your research, please cite:

> Diendorfer AB, Khamina K, Pultar M, Hackl M. miND (miRNA NGS Discovery pipeline): a small RNA-seq analysis pipeline and report generator for microRNA biomarker discovery studies. *F1000Research*. 2022;11:233. doi: [10.12688/f1000research.107833.1](https://doi.org/10.12688/f1000research.107833.1)

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata.

## License

This project is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.
