# Plant Plasmidome Analysis Pipeline

A comprehensive Nextflow pipeline for the large-scale analysis of plant microbiome and plasmidome data. This workflow automates the process from raw sequencing reads to functional annotation and phylogenetic analysis.

## Overview

The pipeline is designed to handle both raw sequencing data and pre-assembled/labeled genomes. It integrates multiple state-of-the-art tools to perform:

1.  **Assembly & QC**: Spades assembly and quality filtering (Raw mode).
2.  **Classification**: Plasmid/virus detection using GeNomad.
3.  **Quality Assessment**: Genome quality check using CheckM.
4.  **Taxonomy**: Taxonomic classification using GTDB-Tk.
5.  **Annotation**:
    *   **Prokka**: General genomic annotation.
    *   **Prodigal**: Gene prediction.
    *   **EggNOG**: Functional annotation (optional).
6.  **Plasmid Analysis**: MOB-suite for plasmid reconstruction and typing.
7.  **Defense Systems**: Comprehensive detection of defense systems using:
    *   **PADLOC**
    *   **DefenseFinder**
    *   **CCTyper** (CRISPR-Cas)
    *   **AMRFinder** (Antimicrobial resistance)
8.  **Phylogeny**: Phylogenetic tree construction using IQ-TREE.

## Prerequisites

*   **Nextflow** (>=23.04.0)
*   **Conda** or **Mamba** (for environment management)
*   **Databases**:
    *   GTDB-Tk database
    *   EggNOG database (if running EggNOG)
    *   Genomad database
    *   CheckM data

## Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/yourusername/Plasmidome_In_Plants.git
    cd Plasmidome_In_Plants
    ```

2.  Ensure all necessary Conda environments are created as specified in `nextflow.config`.

## Usage

The pipeline supports two main input modes: **Raw** and **Labeled**.

### 1. Raw Mode (Start from FASTQ)

Use this mode if you have raw paired-end FASTQ sequencing data.

```bash
nextflow run main.nf \
  --input_mode raw \
  --raw_reads "path/to/Raw_Reads" \
  --outdir "Results"
```

### 2. Labeled Mode (Start from FASTA)

Use this mode if you already have labeled FASTA files (e.g., from a previous run or external source) and a CheckM QA file.

```bash
nextflow run main.nf \
  --input_mode labeled \
  --labeled_fastas "path/to/labeled_fastas" \
  --existing_qa "path/to/qa.tsv" \
  --outdir "Results"
```

For detailed instructions on Labeled Mode, see [LABELED_MODE_USAGE.md](LABELED_MODE_USAGE.md).

### 3. Defense Systems Analysis

The pipeline automatically runs PADLOC, DefenseFinder, and CCTyper. To merge these results into a unified report:

```bash
nextflow run main.nf \
  ... \
  --run_defense_merger true \
  --defense_mapping_file "data/Defense_Systems_Name_List.xlsx"
```

For details on the merger module, see [DEFENSE_MERGER_USAGE.md](DEFENSE_MERGER_USAGE.md).

## Configuration

Key parameters can be configured in `nextflow.config` or passed via command line:

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--length_threshold` | 1000 | Minimum contig length |
| `--coverage_threshold` | 20 | Minimum coverage |
| `--completeness_threshold` | 90.0 | Minimum CheckM completeness |
| `--contamination_threshold` | 5.0 | Maximum CheckM contamination |
| `--run_gtdbtk` | false | Enable GTDB-Tk classification |
| `--run_eggnog` | false | Enable EggNOG annotation |
| `--run_grouped_trees` | false | Build separate trees for sample groups |

## Output Structure

Results are organized in the output directory (default: `Results/`):

```
Results/
├── hq_genomes/           # High-quality genome lists and summaries
├── gtdbtk/               # Taxonomic classification results
├── mobsuite_recon/       # Plasmid reconstruction results
├── mobsuite_typer/       # Plasmid typing results
├── prokka/               # Genome annotations (GFF, FAA, FNA)
├── prodigal/             # Predicted genes and proteins
├── padloc/               # Defense system annotations (PADLOC)
├── defensefinder/        # Defense system annotations (DefenseFinder)
├── cctyper/              # CRISPR-Cas annotations
├── amrfinder/            # AMR gene annotations
├── defense_merger/       # Unified defense system reports
├── iqtree/               # Phylogenetic trees
└── pipeline_info/        # Execution reports (timeline, trace, DAG)
```

## Project Structure

```
├── Collect/          # Data collection scripts
├── Reference/        # Reference documentation
├── Report/           # Analysis reports
├── Result/           # Processed results
├── bin/              # Helper scripts
├── data/             # Data files (e.g., mapping lists)
├── modules/          # Nextflow DSL2 modules
├── main.nf           # Main pipeline script
├── nextflow.config   # Pipeline configuration
└── README.md         # This file
```

## License

MIT License - see LICENSE file for details.
