# Script workflow and cleanup notes

This file documents the current root-level scripts after renaming and light
deduplication. The goal is to keep the original analysis behavior intact while
making the script order and responsibilities easier to follow.

## Naming convention

Root-level analysis scripts use numeric prefixes:

- `001`: raw data download or raw input preparation.
- `101-102`: canonical master-table construction and small table helpers.
- `201-211`: master-table downstream summaries, figures, and tree/iTOL inputs.
- `212-218`: gene-type expansion, PCoA calculation, PCoA plotting, and related diagnostics.
- `219-224`: network table preparation and transfer/plasmid overlap summaries.
- `225-232`: network table cleanup and network plotting.

Script names follow this pattern:

`number_action_object_by_or_with_context.extension`

The action word is intentionally limited:

- `build`: build a core derived table or network from upstream inputs.
- `extract`: extract a smaller table from a larger source.
- `filter`: keep a defined subset of samples, genera, or mappings.
- `calculate`: calculate statistics, matrices, richness, or intersections.
- `analyze`: run an interpretation or association analysis.
- `prepare`: generate files for a downstream tool or later step.
- `merge`: merge records or tables.
- `annotate`: add metadata annotations.
- `plot`: generate figures.

## Current script inventory

| Script | Main role | Default input/output area |
| --- | --- | --- |
| `001_download_ena_PRJEB11584_leaf_http.sh` | Download raw ENA reads | Current working directory |
| `101_build_master_table.py` | Build canonical master table from normalized inputs | `Collect/NCBI_4395_Batch/Master_Table` |
| `102_extract_first_sample_per_species.py` | Extract the first sample for each species | `Collect/.../Master_Table/final` to `Result/.../05_Tree` |
| `201_plot_top_gene_subtypes_by_host_and_contig_type.R` | Top gene subtype bar plots by host and contig type | `Result/.../01_Gene_Distribution_Bar` |
| `202_plot_gene_subtype_stacked_bars_by_host_and_contig_type.R` | Stacked gene-subtype bar plots by host and contig type | `Result/.../02_Gene_Distribution_Stackbar` |
| `203_calculate_gene_richness_by_sample.R` | Richness summaries and alpha-diversity plots | `Result/.../03_Gene_Diversity/Alpha_Diversity` |
| `204_plot_gene_overlap_by_crop_and_contig_type.R` | Venn/Upset-style crop and contig overlap plots | `Result/.../03_Gene_Diversity` |
| `205_plot_taxonomic_gene_heatmaps.R` | Taxonomy-level gene heatmaps | `Result/.../04_Gene_Taxonomy` |
| `206_calculate_genus_level_statistics.py` | Genus-level summary tables for iTOL | `Result/.../05_Tree/Genus_Level/Tree_annotation_file_prepare` |
| `207_prepare_itol_annotation_files.py` | Convert genus statistics to iTOL annotation files | `Result/.../05_Tree/Genus_Level/Tree_visualization` |
| `208_extract_taxonomy_matched_samples.py` | Extract GTDB metadata rows matching target genera | `Result/.../05_Tree/Genus_Level` |
| `209_filter_samples_by_genus.py` | Filter master table rows by selected genera | `Result/.../05_Tree/Genus_Level` |
| `210_filter_sample_genus_mapping.py` | Filter sample-genus map by genus intersection | `Result/.../05_Tree/Genus_Level_Intersection` |
| `211_calculate_genus_intersection.py` | Compute genus intersection across hosts | `Result/.../05_Tree/Genus_Level_Intersection` |
| `212_expand_gene_type_indicator_columns.py` | Expand Defense/AMR/AntiDefense strings into count columns | `Result/.../06_Cluter` |
| `213_calculate_pcoa_by_contig_type.R` | Calculate PCoA separately for chromosome and plasmid | `Result/.../06_Cluter/PCoA_Data` |
| `214_plot_pcoa_by_contig_type.R` | Plot PCoA outputs from script 213 | `Result/.../06_Cluter/PCoA_Data` |
| `215_calculate_pcoa_chromosome_plasmid_overlay.R` | Calculate PCoA for chromosome+plasmid overlay plots | `Result/.../06_Cluter` |
| `216_plot_pcoa_chromosome_plasmid_overlay.R` | Plot chromosome+plasmid overlay PCoA outputs | `Result/.../06_Cluter` |
| `217_analyze_amr_taxonomy_association.R` | Test AMR/taxonomy associations | `Result/.../06_Cluter/amr_merge` |
| `218_analyze_pcoa_cluster_features.R` | Diagnose PCoA cluster drivers | `Result/.../06_Cluter/amr_merge` |
| `219_annotate_network_nodes.py` | Annotate network node tables with master-table metadata | `Result/.../07_Network` |
| `220_find_isolated_network_nodes.py` | Detect nodes not present in edge files | `Result/.../07_Network` |
| `221_merge_network_clusters_with_master_table.py` | Merge MMseq cluster rename map with canonical master table | `Result/.../07_Network/coocc_network` |
| `222_build_contig_linkage_network.py` | Build contig-level linkage network tables | `Result/.../07_Network/transfer_network` |
| `223_prepare_transfer_node_rename_mapping.py` | Match transfer-network nodes to rename map and extract original IDs | `Result/.../07_Network/transfer_network` |
| `224_plot_transferred_plasmid_gene_overlap.py` | Plot transferred vs non-transferred plasmid gene-set overlap | `Result/.../07_Network/transferred_plasmid_gene_overlap` |
| `225_merge_isolate_duplicate_nodes.R` | Merge duplicate isolate-network nodes | `Result/.../07_Network/isolate_network` |
| `226_remove_parallel_network_edges.R` | Remove duplicate undirected isolate-network edges | `Result/.../07_Network/isolate_network` |
| `227_plot_cooccurrence_network.R` | Plot co-occurrence network | `Result/.../07_Network/coocc_network` |
| `228_plot_labeled_cooccurrence_network.R` | Plot labeled co-occurrence network | `Result/.../07_Network/coocc_network` |
| `229_plot_isolate_network.R` | Plot isolate network | `Result/.../07_Network/isolate_network` |
| `230_plot_labeled_isolate_network.R` | Plot labeled isolate network | `Result/.../07_Network/isolate_network` |
| `231_plot_compact_isolate_network.R` | Plot compact isolate network variant | `Result/.../07_Network/isolate_network` |
| `232_plot_transfer_network.R` | Plot transfer network | `Result/.../07_Network/transfer_network` |

## Redundancy removed in this pass

- The former standalone original-ID extraction helper was merged into
  `223_prepare_transfer_node_rename_mapping.py` through `--rename-map-only`.
  This keeps the old helper behavior available without keeping a separate
  script that duplicated the final extraction step.

## Remaining merge candidates

These are not merged yet because each still has slightly different output
contracts or plot styling. They are good candidates for the next cleanup pass.

| Candidate scripts | Suggested consolidation | Reason |
| --- | --- | --- |
| `201_plot_top_gene_subtypes_by_host_and_contig_type.R` + `202_plot_gene_subtype_stacked_bars_by_host_and_contig_type.R` | `201_plot_gene_distribution_suite.R` | Same canonical master table input and closely related gene-distribution figures. |
| `213_calculate_pcoa_by_contig_type.R` + `215_calculate_pcoa_chromosome_plasmid_overlay.R` | Shared PCoA calculation helper or one script with `--mode split/overlay` | Both construct gene matrices and run PCoA; the main difference is the sample unit and output folder layout. |
| `214_plot_pcoa_by_contig_type.R` + `216_plot_pcoa_chromosome_plasmid_overlay.R` | One PCoA plotting script with `--mode split/overlay` | Both read PCoA coordinates and variances, then produce PDF plots. |
| `227_plot_cooccurrence_network.R` to `232_plot_transfer_network.R` | One generic network plotting script with presets | These scripts share most of the same igraph/qgraph layout, color, size, and output-device logic. |
| `225_merge_isolate_duplicate_nodes.R` + `226_remove_parallel_network_edges.R` | `225_clean_isolate_network_tables.R` | Both are small isolate-network table cleanup utilities. |
| `208_extract_taxonomy_matched_samples.py` to `211_calculate_genus_intersection.py` | A tree-preparation wrapper script | These are sequential tree/genus input preparation steps, but each still has a clear standalone input/output file. |

## Notes for future cleanup

- The strongest remaining redundancy is in the network plotting scripts
  (`227-232`). A safe next refactor would extract shared plotting functions
  into one helper file first, then replace the current scripts with thin
  preset wrappers.
- The PCoA scripts should be merged only after confirming which existing
  output folders are still used in figures or reports.
- The root scripts now use project-relative default paths, but command-line
  arguments still allow overriding inputs and outputs.
