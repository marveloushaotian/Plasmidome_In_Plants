#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------- Parameters ----------------
params.raw_reads         = "Raw_Reads"                      // Directory containing paired-end FASTQ files
params.outdir            = "Results"                        // Output directory
params.adapters          = "/scratch/miniforge3/envs/anvio-8/share/bbmap/resources/adapters.fa"

// Sample grouping for separate phylogenetic trees
params.sample_groups     = null                  // Path to tab-separated file: Sample_ID	Group_Name (required when run_grouped_trees=true)
params.run_grouped_trees = false                 // Enable grouped tree construction
params.run_overall_tree  = false                 // Enable overall tree construction (all samples together)

// CheckM version and mode selection
params.checkm_version    = "checkm1"             // "checkm1" (default) or "checkm2" (faster, ML-based)
params.checkm_mode       = "single"              // "single" (incremental, recommended) or "batch" (legacy, only for checkm1)
params.checkm2_db        = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/Database/checkM2/CheckM2_database/uniref100.KO.1"  // CheckM2 database path (without .dmnd extension)

params.genomad_db        = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/genomad_db"
params.gtdbtk_db         = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/gtdbtk_db/release226"

// Optional analyses
params.run_gtdbtk        = false                 // Enable GTDB-Tk classification (disabled by default, very time-consuming)

// Quality & QC thresholds (consumed in modules/config)
params.length_threshold        = 1000
params.coverage_threshold      = 20
params.completeness_threshold  = 90.0
params.contamination_threshold = 5.0

// Bootstrap for IQ-TREE
params.bootstrap_reps = 1000
params.alrt_reps      = 1000

// Prokka parameters
params.prokka_kingdom = 'Bacteria'  // or 'Archaea'
params.prokka_rfam = true

// MOB-suite parameters
params.mobsuite_db = null  // Optional: path to MOB-suite database if not using default

// EggNOG parameters
params.eggnog_db = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/eggnog_db"
params.run_eggnog = false  // Enable EggNOG functional annotation (disabled by default, time-consuming)

// ---------------- Include modules ----------------
include { ASSEMBLY            } from './modules/assembly'
include { QUALITY_FILTER      } from './modules/quality_filter'
include { GENOMAD             } from './modules/genomad'
include { LABEL_CONTIGS       } from './modules/label_contigs'
include { PREP_CHECKM_INPUTS  } from './modules/prep_checkm_inputs'
include { CHECKM_BATCH        } from './modules/checkm_batch'
include { CHECKM_SINGLE       } from './modules/checkm_single'
include { CHECKM2_SINGLE      } from './modules/checkm2_single'
include { MERGE_CHECKM_RESULTS } from './modules/merge_checkm_results'
include { FILTER_HQ_GENOMES   } from './modules/filter_hq_genomes'
include { GTDBTK              } from './modules/gtdbtk'
include { GTDBTK_SPLIT; GTDBTK_CLASSIFY_BATCH; GTDBTK_MERGE } from './modules/gtdbtk_batch'
include { PHYLOGENY           } from './modules/phylogeny'
include { PHYLOGENY_GROUPED   } from './modules/phylogeny_grouped'
include { EGGNOG_DIAMOND; EGGNOG_ANNOTATION } from './modules/eggnog'
include { IQTREE              } from './modules/iqtree'
include { IQTREE_GROUPED      } from './modules/iqtree_grouped'
include { MOBSUITE_TYPER } from './modules/mobsuite'
include { PROKKA_SANITIZE; PROKKA_ANNOTATE } from './modules/prokka'
include { PRODIGAL_PREDICT } from './modules/prodigal'
include { PADLOC_ANNOTATE } from './modules/padloc'
include { DEFENSEFINDER_ANNOTATE } from './modules/defensefinder'
include { CCTYPER_ANNOTATE } from './modules/cctyper'
include { AMRFINDER_ANNOTATE } from './modules/amrfinder'
include { COLLECT_GFF_FILES; COLLECT_PADLOC_RESULTS; COLLECT_DEFENSEFINDER_RESULTS; COLLECT_CCTYPER_RESULTS; MERGE_DEFENSE_SYSTEMS; DEFENSE_SUMMARY } from './modules/defense_merger'
include { EXTRACT_CONTIG_MAPPING } from './modules/extract_contig_mapping'

// ---------------- Workflow ----------------
workflow {

    // 0) Build input channels
    reads_ch = Channel
        .fromFilePairs("${params.raw_reads}/*_{1,2}.fastq.gz")
        .map { sid, files -> tuple(sid.toString(), files) }

    adapters_ch = Channel.value( file(params.adapters) )
    genomad_db_ch = Channel.value( params.genomad_db )
    gtdb_db_ch    = Channel.value( params.gtdbtk_db )
    eggnog_db_ch  = Channel.value( params.eggnog_db )

    // 1) Assembly
    //    ASSEMBLY expects: tuple(id, [R1,R2]), path(adapters)
    ASSEMBLY( reads_ch, adapters_ch )
    // ASSEMBLY.out.contigs => tuple(id, path_contigs)

    // 2) Length + coverage filter
    //    QUALITY_FILTER expects: tuple(id, fasta, reads)
    quality_in_ch = ASSEMBLY.out.contigs.join( reads_ch )
    QUALITY_FILTER( quality_in_ch )
    // QUALITY_FILTER.out.filtered_fasta => tuple(id, path_filtered_fasta)

    // 2.5) Filter out empty filtered fastas (samples with no contigs passing QC)
    //      This prevents downstream errors in GENOMAD, LABEL_CONTIGS, and CHECKM
    QUALITY_FILTER.out.filtered_fasta
        .branch { sample_id, fasta ->
            non_empty: fasta.size() > 0
            empty: true
        }
        .set { filtered_results }

    // Log samples with empty filtered fastas to file (silently)
    filtered_results.empty
        .map { sample_id, fasta -> sample_id }
        .collectFile(name: 'empty_samples_skipped.txt', newLine: true, storeDir: params.outdir)

    // Use non-empty samples for downstream analyses
    filtered_fasta_for_downstream = filtered_results.non_empty

    // 3) GeNomad (on filtered contigs)
    //    GENOMAD expects: tuple(id, fasta), val(genomad_db)
    GENOMAD( filtered_fasta_for_downstream, genomad_db_ch )
    // GENOMAD.out.summary => tuple(id, path_genomad_dir)

    // 4) Label contigs using GeNomad results
    //    LABEL_CONTIGS expects: tuple(id, fasta, genomad_dir)
    GENOMAD.out.summary
        .join( filtered_fasta_for_downstream )                     // (id, genomad_dir, fasta)
        .map { id, genomad_dir, fasta -> tuple(id, fasta, genomad_dir) }
        .set { label_in_ch }

    LABEL_CONTIGS( label_in_ch )
    // LABEL_CONTIGS.out.labeled_fasta => tuple(id, labeled_fasta)

    // ========================================
    // 4.5) Extract contig mapping (sample name -> contig name)
    //      This creates a CSV file mapping sample IDs to their contig names
    // ========================================
    
    // Collect all labeled fasta files
    labeled_fasta_list = LABEL_CONTIGS.out.labeled_fasta
        .map { id, fasta -> fasta }
        .collect()
    
    // Extract contig mapping
    EXTRACT_CONTIG_MAPPING( labeled_fasta_list )
    // EXTRACT_CONTIG_MAPPING.out.mapping_table => path("contig_mapping.csv")

    // ========================================
    // 5-7) CheckM Quality Assessment
    //      Version selection: checkm2 (default, faster) or checkm1 (legacy)
    //      Mode selection (checkm1 only): single or batch
    // ========================================
    
    if (params.checkm_version == "checkm2") {
        // === CheckM2 mode (default) ===
        // CheckM2 is faster and more accurate than CheckM1
        // Uses machine learning and only runs in single-sample mode
        
        println "[INFO] Using CheckM2 for quality assessment (faster, ML-based)"
        
        CHECKM2_SINGLE( LABEL_CONTIGS.out.labeled_fasta )
        
        // Collect all QA results
        qa_all_list = CHECKM2_SINGLE.out.qa_all
            .map { id, qa -> qa }
            .collect()
        
        qa_chro_list = CHECKM2_SINGLE.out.qa_chro
            .map { id, qa -> qa }
            .collect()
        
        // Merge results from all samples
        MERGE_CHECKM_RESULTS( qa_all_list, qa_chro_list )
        
        // Set output channels for downstream
        checkm_qa_all = MERGE_CHECKM_RESULTS.out.qa_all
        checkm_qa_chro = MERGE_CHECKM_RESULTS.out.qa_chro
        
    } else if (params.checkm_version == "checkm1") {
        // === CheckM1 mode (legacy) ===
        println "[INFO] Using CheckM1 for quality assessment (legacy mode)"
        
        if (params.checkm_mode == "single") {
            // === Single-sample mode (incremental) ===
            // Each sample runs CheckM independently
            // Advantages: 
            //   - New samples don't trigger re-running old samples
            //   - Better for incremental analysis
            //   - Nextflow cache works per sample
            
            CHECKM_SINGLE( LABEL_CONTIGS.out.labeled_fasta )
            
            // Collect all QA results
            qa_all_list = CHECKM_SINGLE.out.qa_all
                .map { id, qa -> qa }
                .collect()
            
            qa_chro_list = CHECKM_SINGLE.out.qa_chro
                .map { id, qa -> qa }
                .collect()
            
            // Merge results from all samples
            MERGE_CHECKM_RESULTS( qa_all_list, qa_chro_list )
            
            // Set output channels for downstream
            checkm_qa_all = MERGE_CHECKM_RESULTS.out.qa_all
            checkm_qa_chro = MERGE_CHECKM_RESULTS.out.qa_chro
            
        } else {
            // === Batch mode (legacy) ===
            // All samples run CheckM together
            // Advantages:
            //   - Slightly faster for small datasets
            //   - Single CheckM invocation
            
            labeled_fasta_list = LABEL_CONTIGS.out.labeled_fasta
                .map { id, f -> f }
                .collect()
            
            PREP_CHECKM_INPUTS( labeled_fasta_list )
            
            CHECKM_BATCH(
                PREP_CHECKM_INPUTS.out.bins_all_dir,
                PREP_CHECKM_INPUTS.out.bins_chro_dir
            )
            
            // Set output channels for downstream
            checkm_qa_all = CHECKM_BATCH.out.qa_all
            checkm_qa_chro = CHECKM_BATCH.out.qa_chro
        }
    } else {
        error "ERROR: params.checkm_version must be 'checkm2' or 'checkm1', got: ${params.checkm_version}"
    }
    
    // 8) Build HQ/LQ lists from QA tables
    checkm_qa_all
        .mix( checkm_qa_chro )
        .collect()
        .set { all_qa_files_ch }

    FILTER_HQ_GENOMES( all_qa_files_ch )
    // FILTER_HQ_GENOMES.out.hq_list => path("hq_genomes.txt")

    // ========================================
    // 9) Create unified HQ genomes dataset
    //    This is the SINGLE SOURCE for all HQ genome files
    //    All downstream analyses should reference this channel
    // ========================================
    hq_ids_ch = FILTER_HQ_GENOMES.out.hq_list
        .splitText()
        .map { it.trim() }
        .filter { it }

    // Create HQ genomes channel based on CheckM mode
    if (params.checkm_mode == "single") {
        // In single mode, use labeled_fasta files directly
        hq_genomes_ch = LABEL_CONTIGS.out.labeled_fasta
            .join( hq_ids_ch.map { id -> tuple(id, id) } )  // Keep only HQ samples
            .map { id, fasta, hq_id -> tuple(id, fasta) }
    } else {
        // In batch mode, use PREP_CHECKM_INPUTS output
        hq_genomes_ch = PREP_CHECKM_INPUTS.out.bins_all_dir
            .combine( hq_ids_ch )                                          // (dir, id)
            .map { dir, id -> tuple(id, file("${dir}/${id}.fasta")) }      // tuple(id, fasta)
    }
    // hq_genomes_ch is now the unified channel: tuple(genome_id, fasta_file)

    // ========================================
    // 10) GTDB-Tk classification on HQ genomes (optional, disabled by default)
    //     Requires: list of FASTA files
    //     Enable with: --run_gtdbtk true
    // ========================================
    if (params.run_gtdbtk) {
        println "[INFO] GTDB-Tk classification enabled (batch mode: 50 genomes/batch)"

        hq_fasta_list = hq_genomes_ch
            .map { id, fasta -> fasta }
            .collect()

        // Batch processing: split -> classify in parallel -> merge
        GTDBTK_SPLIT( hq_fasta_list )
        GTDBTK_CLASSIFY_BATCH( 
            GTDBTK_SPLIT.out.batches.flatten(),
            gtdb_db_ch 
        )
        GTDBTK_MERGE( 
            GTDBTK_CLASSIFY_BATCH.out.summary.collect(),
            GTDBTK_CLASSIFY_BATCH.out.alignments.collect()
        )
        
        // Output channels (compatible with downstream processes)
        gtdbtk_summary = GTDBTK_MERGE.out.summary
        gtdbtk_alignments = GTDBTK_MERGE.out.alignments
    } else {
        println "[INFO] GTDB-Tk classification disabled (set --run_gtdbtk true to enable)"
    }

    // ========================================
    // 11) MOB-suite analysis for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    // MOBSUITE_RECON( hq_genomes_ch )  // Disabled
    MOBSUITE_TYPER( hq_genomes_ch )

    // ========================================
    // 12) Prokka annotation for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    PROKKA_SANITIZE( hq_genomes_ch )
    PROKKA_ANNOTATE( PROKKA_SANITIZE.out.sanitized_fasta )

    // ========================================
    // 13) Prodigal gene prediction for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    //     This is the SOURCE for protein/gene sequences
    // ========================================
    PRODIGAL_PREDICT( hq_genomes_ch )

    // ========================================
    // 14) DefenseFinder annotation for HQ genomes
    //     Requires: tuple(id, proteins.faa) from Prodigal
    // ========================================
    DEFENSEFINDER_ANNOTATE( PRODIGAL_PREDICT.out.proteins )

    // ========================================
    // 15) CCTyper CRISPR-Cas annotation for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    CCTYPER_ANNOTATE( hq_genomes_ch )

    // ========================================
    // 16) AMRFinder resistance gene annotation for HQ genomes
    //     Requires: tuple(id, proteins.faa) from Prodigal
    // ========================================
    AMRFINDER_ANNOTATE( PRODIGAL_PREDICT.out.proteins )

    // ========================================
    // 17) PADLOC defense systems annotation for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    PADLOC_ANNOTATE( hq_genomes_ch )

    // ========================================
    // 18-19) Phylogenetic tree construction (requires GTDB-Tk)
    //        Only runs if params.run_gtdbtk = true
    // ========================================

    if (params.run_gtdbtk) {
        // Tree construction modes:
        // 1. Overall tree: all samples together (params.run_overall_tree = true)
        // 2. Grouped trees: separate tree for each group (params.run_grouped_trees = true)
        // 3. Both: run both overall and grouped trees

        if (params.run_overall_tree) {
            // Original single tree construction for all samples
            PHYLOGENY( gtdbtk_alignments, FILTER_HQ_GENOMES.out.hq_list )
            IQTREE( PHYLOGENY.out.filtered_alignments )
        }

        if (params.run_grouped_trees) {
            // Validate that sample_groups file is provided
            if (params.sample_groups == null) {
                error "ERROR: params.sample_groups must be specified when params.run_grouped_trees = true"
            }

            // Load sample grouping information only when needed
            // Creates tuples: (sample_id, group_name)
            sample_groups_pairs = Channel.fromPath( params.sample_groups, checkIfExists: true )
                .splitCsv( header: false, sep: '\t', comment: '#' )
                .filter { row -> row.size() >= 2 && row[0].trim() && row[1].trim() }  // Filter empty lines
                .map { row -> tuple(row[0].trim(), row[1].trim()) }  // (sample_id, group_name)

            // Create channel combining HQ genomes with their group information
            // Join HQ genome IDs with their groups, then group by group_name
            hq_genomes_with_groups_ch = hq_ids_ch
                .combine( sample_groups_pairs )  // (genome_id, sample_id, group_name)
                .filter { genome_id, sample_id, group_name -> genome_id == sample_id }  // Keep only matching IDs
                .map { genome_id, sample_id, group_name -> tuple(group_name, genome_id) }  // (group_name, genome_id)
                .groupTuple()  // Group genome IDs by group name: (group_name, [genome_id1, genome_id2, ...])

            // Grouped tree construction: separate trees for each group
            PHYLOGENY_GROUPED(
                gtdbtk_alignments,
                FILTER_HQ_GENOMES.out.hq_list,
                hq_genomes_with_groups_ch
            )

            // IQ-TREE for each group
            IQTREE_GROUPED( PHYLOGENY_GROUPED.out.filtered_alignments )
        }
    }

    // ========================================
    // 20) EggNOG functional annotation for HQ genomes (optional)
    //     Two-step process: 1) Diamond search 2) Annotation
    // ========================================
    if (params.run_eggnog) {
        println "[INFO] EggNOG functional annotation enabled"
        
        // Collect all prokka .faa files
        prokka_faa_files = PROKKA_ANNOTATE.out.proteins
            .map { id, faa -> faa }
            .collect()
        
        // Step 1: Diamond search
        EGGNOG_DIAMOND( prokka_faa_files, eggnog_db_ch )
        
        // Step 2: Annotation
        EGGNOG_ANNOTATION( EGGNOG_DIAMOND.out.diamond_dir, eggnog_db_ch )
    } else {
        println "[INFO] EggNOG annotation disabled (set --run_eggnog true to enable)"
    }

    // ========================================
    // 21) Merge defense systems annotations (always enabled)
    //     Combines results from PADLOC, DefenseFinder, and CCTyper
    //     Mapping file: Defense_Systems_Name_List.xlsx (in project root)
    // ========================================
    println "[INFO] Merging defense systems annotations from PADLOC, DefenseFinder, and CCTyper"

    // Collect GFF files from Prodigal (defense tools use Prodigal annotations)
    prodigal_gff_files = PRODIGAL_PREDICT.out.gff
        .map { id, gff -> gff }
        .collect()

    COLLECT_GFF_FILES( prodigal_gff_files )

    // Collect PADLOC results
    padloc_results = PADLOC_ANNOTATE.out.results
        .map { id, tsv -> tsv }
        .collect()

    COLLECT_PADLOC_RESULTS( padloc_results )

    // Collect DefenseFinder results
    defensefinder_results = DEFENSEFINDER_ANNOTATE.out.results
        .map { id, tsv -> tsv }
        .collect()

    COLLECT_DEFENSEFINDER_RESULTS( defensefinder_results )

    // Collect CCTyper results (cas_operons.tab files)
    cctyper_results = CCTYPER_ANNOTATE.out.cas_operons
        .map { id, tsv -> tsv }
        .collect()

    COLLECT_CCTYPER_RESULTS( cctyper_results )

    // Load mapping file (fixed path in project root)
    mapping_file_ch = Channel.fromPath("Defense Systems Name List.xlsx", checkIfExists: true)

    // Merge all defense systems annotations
    MERGE_DEFENSE_SYSTEMS(
        COLLECT_GFF_FILES.out.merged_gff,
        COLLECT_PADLOC_RESULTS.out.merged_padloc,
        COLLECT_DEFENSEFINDER_RESULTS.out.merged_defensefinder,
        COLLECT_CCTYPER_RESULTS.out.merged_cctyper,
        mapping_file_ch
    )

    // Generate summary report
    DEFENSE_SUMMARY( MERGE_DEFENSE_SYSTEMS.out.defense_summary )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Status: ${ workflow.success ? 'OK' : 'FAILED' }"
    println "Duration: $workflow.duration"
}