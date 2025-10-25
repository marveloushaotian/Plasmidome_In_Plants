#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------- Parameters ----------------
params.project_id        = "TEST"                           // Project identifier for output organization
params.raw_reads         = "Raw_Reads"                      // Directory containing paired-end FASTQ files
params.outdir            = "Results/${params.project_id}"   // Output directory
params.adapters          = "/scratch/miniforge3/envs/anvio-8/share/bbmap/resources/adapters.fa"

// Sample grouping for separate phylogenetic trees
params.sample_groups     = null                  // Path to tab-separated file: Sample_ID	Group_Name (required when run_grouped_trees=true)
params.run_grouped_trees = false                 // Enable grouped tree construction
params.run_overall_tree  = false                 // Enable overall tree construction (all samples together)

// CheckM mode selection
params.checkm_mode       = "single"              // "single" (incremental, recommended) or "batch" (legacy)

params.genomad_db        = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/genomad_db"
params.gtdbtk_db         = "/projects/somicrobiology/data/CRBC_CropMicrobiomeBacteria/ketao/gtdbtk_db/release226"

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
params.eggnog_db = "/dev/shm/eggnog_db/"
params.run_eggnog = true                    // Enable EggNOG functional annotation
params.eggnog_protein_source = "prodigal"   // Options: "prodigal" (default, faster) or "prokka" (more detailed) 

// ---------------- Include modules ----------------
include { ASSEMBLY            } from './modules/assembly'
include { QUALITY_FILTER      } from './modules/quality_filter'
include { GENOMAD             } from './modules/genomad'
include { LABEL_CONTIGS       } from './modules/label_contigs'
include { PREP_CHECKM_INPUTS  } from './modules/prep_checkm_inputs'
include { CHECKM_BATCH        } from './modules/checkm_batch'
include { CHECKM_SINGLE       } from './modules/checkm_single'
include { MERGE_CHECKM_RESULTS } from './modules/merge_checkm_results'
include { FILTER_HQ_GENOMES   } from './modules/filter_hq_genomes'
include { GTDBTK              } from './modules/gtdbtk'
include { PHYLOGENY           } from './modules/phylogeny'
include { PHYLOGENY_GROUPED   } from './modules/phylogeny_grouped'
include { EGGNOG_DIAMOND; EGGNOG_ANNOTATION } from './modules/eggnog'
include { IQTREE              } from './modules/iqtree'
include { IQTREE_GROUPED      } from './modules/iqtree_grouped'
include { MOBSUITE_RECON; MOBSUITE_TYPER } from './modules/mobsuite'
include { PROKKA_SANITIZE; PROKKA_ANNOTATE } from './modules/prokka'
include { PRODIGAL_PREDICT } from './modules/prodigal'
include { PADLOC_ANNOTATE } from './modules/padloc'
include { DEFENSEFINDER_ANNOTATE } from './modules/defensefinder'
include { CCTYPER_ANNOTATE } from './modules/cctyper'
include { AMRFINDER_ANNOTATE } from './modules/amrfinder'

// ---------------- Workflow ----------------
workflow {

    // 0) Build input channels
    reads_ch = Channel
        .fromFilePairs("${params.raw_reads}/*_{1,2}.fastq.gz")
        .map { sid, files -> tuple(sid.toString(), files) }

    eggnog_db_ch = Channel.value( params.eggnog_db )
    adapters_ch = Channel.value( file(params.adapters) )
    genomad_db_ch = Channel.value( params.genomad_db )
    gtdb_db_ch    = Channel.value( params.gtdbtk_db )

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

    // Log samples with empty filtered fastas
    filtered_results.empty
        .subscribe { sample_id, fasta ->
            log.warn "⚠️  Sample ${sample_id}: filtered fasta is empty (no contigs passed QC thresholds: length>=${params.length_threshold}, coverage>=${params.coverage_threshold}). Skipping all downstream analyses for this sample."
        }

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
    // 5-7) CheckM Quality Assessment
    //      Two modes available:
    //      - single: Run CheckM per sample (incremental, recommended for large-scale)
    //      - batch:  Run CheckM in batch (legacy, faster for small datasets)
    // ========================================
    
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
    
    // 8) Build HQ/LQ lists from QA tables
    checkm_qa_all
        .mix( checkm_qa_chro )
        .collect()
        .set { all_qa_files_ch }

    FILTER_HQ_GENOMES( all_qa_files_ch )
    // FILTER_HQ_GENOMES.out.hq_list => path("hq_genomes.txt")

    // ========================================
    // 8) Create unified HQ genomes dataset
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
    // 9) GTDB-Tk classification on HQ genomes (batch)
    //    Requires: list of FASTA files
    // ========================================
    hq_fasta_list = hq_genomes_ch
        .map { id, fasta -> fasta }
        .collect()

    GTDBTK( hq_fasta_list, gtdb_db_ch )
    // GTDBTK.out.alignments => path("gtdbtk_output/align")

    // ========================================
    // 10) MOB-suite analysis for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    MOBSUITE_RECON( hq_genomes_ch )
    MOBSUITE_TYPER( hq_genomes_ch )

    // ========================================
    // 11) Prokka annotation for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    PROKKA_SANITIZE( hq_genomes_ch )
    PROKKA_ANNOTATE( PROKKA_SANITIZE.out.sanitized_fasta )

    // ========================================
    // 12) Prodigal gene prediction for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    //     This is the SOURCE for protein/gene sequences
    // ========================================
    PRODIGAL_PREDICT( hq_genomes_ch )

    // ========================================
    // 13) DefenseFinder annotation for HQ genomes
    //     Requires: tuple(id, proteins.faa) from Prodigal
    // ========================================
    DEFENSEFINDER_ANNOTATE( PRODIGAL_PREDICT.out.proteins )

    // ========================================
    // 14) CCTyper CRISPR-Cas annotation for HQ genomes
    //     Requires: tuple(id, genes.fna) from Prodigal
    // ========================================
    CCTYPER_ANNOTATE( PRODIGAL_PREDICT.out.genes )

    // ========================================
    // 15) AMRFinder resistance gene annotation for HQ genomes
    //     Requires: tuple(id, proteins.faa) from Prodigal
    // ========================================
    AMRFINDER_ANNOTATE( PRODIGAL_PREDICT.out.proteins )

    // ========================================
    // 16) PADLOC defense systems annotation for HQ genomes
    //     Requires: tuple(id, fasta) - uses unified hq_genomes_ch
    // ========================================
    PADLOC_ANNOTATE( hq_genomes_ch )

    // ========================================
    // 17) Filter alignments by HQ ids (for tree building)
    //     Requires: GTDBTK alignments + HQ list + Group info
    // ========================================

    // ========================================
    // Tree construction modes:
    // 1. Overall tree: all samples together (params.run_overall_tree = true)
    // 2. Grouped trees: separate tree for each group (params.run_grouped_trees = true)
    // 3. Both: run both overall and grouped trees
    // ========================================
    
    if (params.run_overall_tree) {
        // Original single tree construction for all samples
        PHYLOGENY( GTDBTK.out.alignments, FILTER_HQ_GENOMES.out.hq_list )
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
            GTDBTK.out.alignments,
            FILTER_HQ_GENOMES.out.hq_list,
            hq_genomes_with_groups_ch
        )

        // IQ-TREE for each group
        IQTREE_GROUPED( PHYLOGENY_GROUPED.out.filtered_alignments )
    }

    // ========================================
    // 19) EggNOG functional annotation for HQ genomes (optional)
    //     Protein source can be configured via params.eggnog_protein_source
    //     - "prodigal": Use Prodigal proteins (default, faster, always available)
    //     - "prokka": Use Prokka proteins (more detailed gene names)
    // ========================================
    if (params.run_eggnog) {
        // Select protein source based on parameter
        if (params.eggnog_protein_source == "prokka") {
            // Use Prokka annotated proteins
            eggnog_proteins_ch = PROKKA_ANNOTATE.out.proteins
                .map { id, f -> tuple(id, f) }
            
            println "[INFO] EggNOG will use Prokka proteins (detailed gene names)"
        } else {
            // Use Prodigal predicted proteins (default)
            eggnog_proteins_ch = PRODIGAL_PREDICT.out.proteins
                .map { id, f -> tuple(id, f) }
            
            println "[INFO] EggNOG will use Prodigal proteins (faster, default)"
        }
        
        // Run eggNOG diamond search
        EGGNOG_DIAMOND( eggnog_proteins_ch, eggnog_db_ch )
        
        // Run eggNOG annotation
        EGGNOG_ANNOTATION( EGGNOG_DIAMOND.out.diamond_results, eggnog_db_ch )
    } else {
        println "[INFO] EggNOG annotation is disabled (set params.run_eggnog=true to enable)"
    }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Status: ${ workflow.success ? 'OK' : 'FAILED' }"
    println "Duration: $workflow.duration"
}
