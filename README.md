# Plasmidome In Plants

This repository contains the codebase for our plant plasmidome research project. The project focuses on processing raw sequencing data from plant-associated microbial communities and following the analysis through genome reconstruction, plasmidome-related annotation, defense system characterization, resistance profiling, taxonomy, and downstream statistical analysis and visualization.

In practice, this repository covers two broad stages of the project:

- primary processing from raw sequencing reads to annotated genome-level result tables
- downstream table-based analysis, integration, statistics, visualization, and network analysis

## Project Scope

The project is designed to study plant-associated bacterial genomes and plasmid-related features from raw paired-end sequencing data. The workflow starts with raw read processing and continues through assembly, contig classification, genome quality control, annotation, plasmid typing, defense system identification, and generation of summary tables for downstream biological interpretation.

This means the repository is not just a pipeline repository. It is the main code repository for the full project, including both the reproducible workflow layer and the scripts used for follow-up biological analysis.

## Main Workflow

The core end-to-end workflow is implemented in the Nextflow directories under [`workflows/nextflow`](workflows/nextflow). As a representative example, [`workflows/nextflow/NEXTFLOW`](workflows/nextflow/NEXTFLOW) shows the canonical structure of the pipeline.

At a high level, the Nextflow workflow does the following:

1. reads paired-end raw FASTQ files
2. performs assembly and quality filtering
3. uses GeNomad to label contigs
4. runs CheckM or CheckM2 for genome quality assessment
5. extracts high-quality genomes for downstream analysis
6. performs genome annotation with Prokka and Prodigal
7. types plasmid-related sequences with MOB-suite
8. detects defense and resistance features with PADLOC, DefenseFinder, CCTyper, and AMRFinder
9. optionally performs GTDB-Tk taxonomic classification, phylogenetic analysis, and EggNOG functional annotation
10. merges defense system results into final summary tables

So the pipeline covers the project from raw sequencing input all the way to structured annotation outputs and summary tables.

## Why There Are Multiple Nextflow Directories

This repository contains four related workflow directories, all grouped under [`workflows/nextflow`](workflows/nextflow):

- [`workflows/nextflow/NEXTFLOW`](workflows/nextflow/NEXTFLOW)
- [`workflows/nextflow/NEXTFLOW_PIPELINE`](workflows/nextflow/NEXTFLOW_PIPELINE)
- [`workflows/nextflow/NEXTFLOW_Computerome`](workflows/nextflow/NEXTFLOW_Computerome)
- [`workflows/nextflow/NEXFFLOW_HALF`](workflows/nextflow/NEXFFLOW_HALF)

These directories reflect the history of the project as the workflow was migrated and adapted across different compute environments and servers during the course of analysis. They should be understood as environment-specific or stage-specific workflow variants rather than independent scientific projects.

For explaining the workflow structure, `NEXTFLOW` is the clearest baseline example. In actual use, a different directory may be preferred depending on the target server and environment configuration.

## Downstream Analysis Scripts

After the Nextflow workflow generates annotation tables and intermediate result files, the master table is assembled with [`101_build_master_table.py`](101_build_master_table.py). The canonical input and output files for this table are kept under [`Collect/NCBI_4395_Batch/Master_Table`](Collect/NCBI_4395_Batch/Master_Table), with the frequently used final table at [`Collect/NCBI_4395_Batch/Master_Table/final/05_master_contig_annotation_table.csv`](Collect/NCBI_4395_Batch/Master_Table/final/05_master_contig_annotation_table.csv).

The numbered Python and R scripts in the repository root are then used for downstream analysis.

These scripts are used for tasks such as:

- calculating diversity, richness, and subtype composition summaries
- preparing taxonomic summaries and comparative matrices
- generating PCoA, heatmaps, stacked bar plots, Venn-style comparisons, and other figures
- preparing iTOL annotations and phylogeny-linked summaries
- constructing and plotting co-occurrence, isolate, and transfer networks

In other words, the numbered scripts are the project’s analysis layer for turning workflow outputs into biological summaries, comparative analyses, and figures.

## Repository Structure

The repository is organized around the full project workflow rather than a single software package.

### Workflow directories

- [`workflows/nextflow/NEXTFLOW`](workflows/nextflow/NEXTFLOW): canonical example of the main pipeline layout
- [`workflows/nextflow/NEXTFLOW_PIPELINE`](workflows/nextflow/NEXTFLOW_PIPELINE): alternate workflow organization used during development
- [`workflows/nextflow/NEXTFLOW_Computerome`](workflows/nextflow/NEXTFLOW_Computerome): Computerome-adapted workflow version
- [`workflows/nextflow/NEXFFLOW_HALF`](workflows/nextflow/NEXFFLOW_HALF): partial or intermediate workflow snapshot

### Root-level analysis scripts

- `001_*.sh`: data download and raw input preparation
- [`101_build_master_table.py`](101_build_master_table.py): canonical master table builder from the four normalized input files to the final `05_master_contig_annotation_table.csv`
- [`102_extract_first_sample_per_species.py`](102_extract_first_sample_per_species.py): helper extraction script retained from the early table-analysis layer
- `201_*.R` onward: statistical analysis, comparative summaries, ordination, visualization, and figure generation
- `217_*.py` onward plus `226_*.R` to `232_*.R`: graph and network annotation, network construction, and network plotting

### Project-linked directories

- [`Collect`](Collect): collected inputs and project materials
- [`Reference`](Reference): reference files and supporting materials
- [`Report`](Report): report outputs
- [`Result`](Result): project results and exports

Depending on the machine, these may be symlinks to external project storage.

## How To Use The Workflow

For someone new to the project, the recommended way to understand and use the repository is:

1. identify which directory under `workflows/nextflow` matches the target compute environment
2. inspect `main.nf`, `nextflow.config`, `modules/`, and `bin/` in that workflow directory
3. prepare raw paired-end reads in the expected naming format
4. confirm environment paths, databases, and required resource files for the selected workflow directory
5. run the Nextflow workflow to generate annotation and summary tables
6. use the numbered Python and R scripts in the repository root for downstream table processing, statistical analysis, and figure generation

For understanding the logic of the pipeline itself, start with:

- [`workflows/nextflow/NEXTFLOW/main.nf`](workflows/nextflow/NEXTFLOW/main.nf)
- [`workflows/nextflow/NEXTFLOW/nextflow.config`](workflows/nextflow/NEXTFLOW/nextflow.config)
- [`workflows/nextflow/NEXTFLOW/modules`](workflows/nextflow/NEXTFLOW/modules)

## Contribution And Project Roles

This project is led by Tao Ke.

The workflow development, adaptation, and downstream analyses in this repository were supported by Haozhe and Pau. We thank all contributors to the project for pipeline development, analysis support, troubleshooting, and scientific discussion.

## Contact

Correspondence:

- [Contact name / email to be added]

## License

MIT License. See [`LICENSE`](LICENSE).
