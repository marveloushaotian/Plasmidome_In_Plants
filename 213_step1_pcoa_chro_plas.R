#!/usr/bin/env Rscript

# =============================================================================
# Step 1 (Chromosome+Plasmid overlay): PCoA Calculation by Host and Contig Type
# Description:
#   - From expanded gene type profiles, compute PCoA for ONE gene type
#     (Defense, AMR, or AntiDefense) per run.
#   - Rows in the PCoA matrix are "Sample_ID × Contig_Type2" combinations.
#   - Output per-gene-type PCoA coordinates and variance explained for use
#     by step-2 plotting with chromosome and plasmid in the same panel.
#
# Usage:
#   Rscript 213_step1_pcoa_chro_plas.R \
#       -i <input.csv> \
#       -t <type> \
#       -o <output_prefix> \
#       [-c]
#
# Arguments:
#   -i / --input      : Input CSV file path (expanded gene type table) for
#                       one gene type (defense, amr, or antidefense).
#   -t / --type       : Gene type: 'defense', 'amr', or 'antidefense'.
#   -o / --output     : Output file prefix (default: PCoA_chro_plas).
#   -c / --use-counts : Use count numbers (Bray-Curtis distance). If not
#                       set, use binary presence/absence (Jaccard distance).
#
# Notes:
#   - No PERMANOVA is performed here.
#   - Coordinates and variance are written to input-dir subfolders:
#       plas_chro_coordinates and plas_chro_variance
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(ape)
  library(argparse)
})

# ------------------------- 1. Parse arguments ------------------------------- #

parser <- ArgumentParser(
  description = "PCoA calculation for chromosome+plasmid overlay design."
)
parser$add_argument(
  "-i", "--input",
  required = TRUE,
  help = "Input CSV file path (expanded gene type table) for one gene type."
)
parser$add_argument(
  "-t", "--type",
  required = TRUE,
  choices = c("defense", "amr", "antidefense"),
  help = "Gene type: 'defense', 'amr', or 'antidefense'."
)
parser$add_argument(
  "-o", "--output",
  default = "PCoA_chro_plas",
  help = "Output prefix for chromosome+plasmid overlay design (default: PCoA_chro_plas)."
)
parser$add_argument(
  "-c", "--use-counts",
  action = "store_true",
  help = "Use count numbers (Bray-Curtis distance). If not set, use binary presence/absence (Jaccard)."
)

args <- parser$parse_args()

# ------------------------- 2. Configuration -------------------------------- #

type_config <- list(
  defense = list(
    allowed_cols = strsplit("Defense_3HP,Defense_6A_MBL,Defense_AVAST type I,Defense_AVAST type II,Defense_AVAST type III,Defense_AVAST type IV,Defense_AVAST type V,Defense_Abi2,Defense_AbiA,Defense_AbiAlpha,Defense_AbiC,Defense_AbiD,Defense_AbiE,Defense_AbiG,Defense_AbiH,Defense_AbiJ,Defense_AbiK,Defense_AbiL,Defense_AbiO/Nhi,Defense_AbiP2,Defense_AbiQ,Defense_AbiR,Defense_AbiU,Defense_AbiV,Defense_AbiZ,Defense_Aditi,Defense_ApeA,Defense_Azaca,Defense_BREX other,Defense_BREX type I,Defense_BREX type II,Defense_BREX type III,Defense_BREX type VI,Defense_Belisama,Defense_Borvo,Defense_Brig1,Defense_BstA,Defense_Butters gp30-gp31,Defense_Butters gp57r,Defense_CARD-NLR-Endonuclease,Defense_CARD-NLR-Phospho,Defense_CARD-NLR-Subtilase,Defense_CARD-NLR-like,Defense_CBASS other,Defense_CBASS type I,Defense_CBASS type II,Defense_CBASS type III,Defense_CBASS type IV,Defense_CRISPR-Cas adaptation,Defense_CRISPR-Cas other,Defense_CRISPR-Cas type I,Defense_CRISPR-Cas type I-G,Defense_CRISPR-Cas type I-III,Defense_CRISPR-Cas type II,Defense_CRISPR-Cas type III,Defense_CRISPR-Cas type IV,Defense_CRISPR-Cas type V,Defense_CRISPR-Cas type VI,Defense_CapRel,Defense_Ceres,Defense_Cernunnos,Defense_Charlie gp32,Defense_CoCoNut I-A,Defense_CoCoNut I-B,Defense_CoCoNut I-C,Defense_CoCoNut II,Defense_CoCoNut III-A,Defense_CoCoNut III-B,Defense_DISARM type I,Defense_DISARM type II,Defense_DMS other,Defense_DRT class I,Defense_DRT class II,Defense_DRT class III,Defense_DRT other,Defense_DRT type 2,Defense_DRT type 3,Defense_DRT type 4,Defense_DRT type 5,Defense_DRT type 6,Defense_DRT type 7,Defense_DRT type 8,Defense_DRT type 9,Defense_DS-1,Defense_DS-10,Defense_DS-11,Defense_DS-12B,Defense_DS-14,Defense_DS-16,Defense_DS-17,Defense_DS-18,Defense_DS-19,Defense_DS-2,Defense_DS-20,Defense_DS-22,Defense_DS-23,Defense_DS-24,Defense_DS-29,Defense_DS-31,Defense_DS-32,Defense_DS-33,Defense_DS-35,Defense_DS-37,Defense_DS-39,Defense_DS-4,Defense_DS-41,Defense_DS-44,Defense_DS-45,Defense_DS-6,Defense_DS-7,Defense_DS-8,Defense_DS-9,Defense_DUF4238,Defense_Damona,Defense_DarTG,Defense_Dazbog,Defense_Detocs,Defense_Divona,Defense_DndABCDE,Defense_DndFGH,Defense_Dodola,Defense_Dpd,Defense_Druantia other,Defense_Druantia type I,Defense_Druantia type II,Defense_Druantia type III,Defense_Druantia type IV,Defense_Dsr1,Defense_Dsr2,Defense_Dynamins,Defense_Eleos,Defense_Epona,Defense_Erebus,Defense_Esos,Defense_GAPS1,Defense_GAPS2,Defense_GAPS4,Defense_GIY-YIG,Defense_Gabija,Defense_Gasdermin,Defense_Geb,Defense_HEC-01,Defense_HEC-02,Defense_HEC-03,Defense_HEC-04,Defense_HEC-05,Defense_HEC-06,Defense_HEC-07,Defense_HEC-08,Defense_HEC-09,Defense_HNH,Defense_HP,Defense_Hachiman,Defense_Hachiman other,Defense_Hachiman type I,Defense_Hachiman type II,Defense_Helicase-DUF2291,Defense_Her-DUF,Defense_Her-SIR,Defense_Hhe,Defense_Hma,Defense_Hna,Defense_Hok/Sok,Defense_Hydrolase-TM,Defense_Hypnos,Defense_ISG15-like,Defense_ISG15-like other,Defense_IetAS,Defense_JukAB,Defense_Kiwa,Defense_Lamassu Amidase,Defense_Lamassu Cap4-nuclease,Defense_Lamassu DS-30,Defense_Lamassu FMO,Defense_Lamassu Hydrolase,Defense_Lamassu Hydrolase-Protease,Defense_Lamassu Lipase,Defense_Lamassu Mrr,Defense_Lamassu PDDEXK,Defense_Lamassu Sir2,Defense_Lamassu family,Defense_Lanthiphage,Defense_Lit,Defense_MADS,Defense_MADS3-4,Defense_MMB gp29-gp30,Defense_MazEF,Defense_Menshen,Defense_Menshen other,Defense_Mokosh type I,Defense_Mokosh type II,Defense_Mza,Defense_Mza other,Defense_NLR_like,Defense_Nantosuelta,Defense_Nemetona,Defense_NixI,Defense_Ogmios,Defense_Old,Defense_Olokun,Defense_Oshun,Defense_PARIS type I,Defense_PARIS type II,Defense_PD-Lambda-1,Defense_PD-Lambda-2,Defense_PD-Lambda-5,Defense_PD-T4-1,Defense_PD-T4-10,Defense_PD-T4-2,Defense_PD-T4-3,Defense_PD-T4-4,Defense_PD-T4-5,Defense_PD-T4-6,Defense_PD-T4-7,Defense_PD-T4-8,Defense_PD-T4-9,Defense_PD-T7-1,Defense_PD-T7-3,Defense_PD-T7-4,Defense_PD-T7-5,Defense_Panchino gp28,Defense_PbeABCD,Defense_PfiAT,Defense_Phrann gp29-gp30,Defense_PifA,Defense_PifAC,Defense_Ppl,Defense_Prithvi,Defense_Prometheus,Defense_PrrC,Defense_Pseudo-CoCoNut,Defense_PsyrTA,Defense_Pycsar,Defense_Pycsar other,Defense_QatABCD,Defense_QatABCD other,Defense_RADAR,Defense_RL,Defense_RM type I,Defense_RM type II,Defense_RM type II-2,Defense_RM type III,Defense_RM type IV,Defense_RM type IV-1,Defense_Retron type I,Defense_Retron type II,Defense_Retron type III,Defense_Retron type IV,Defense_Retron type V,Defense_Retron type VI,Defense_Retron type VII,Defense_Retron type XI,Defense_Retron type XII,Defense_Retron type XIII,Defense_RloC,Defense_RnlAB,Defense_RosmerTA,Defense_SDIC3,Defense_SEFIR,Defense_Septu,Defense_Shango,Defense_Shedu,Defense_ShosTA,Defense_Sirona,Defense_SoFic,Defense_SpbK,Defense_SspABCD,Defense_SspE,Defense_SspFGH,Defense_Stk2,Defense_Sucellos,Defense_TIR-III,Defense_TIR-IV,Defense_TIR-NLR,Defense_Taranis,Defense_TerY-P,Defense_TerY-P other,Defense_TgvAB,Defense_Thoeris other,Defense_Thoeris type I,Defense_Thoeris type II,Defense_Thoeris type III,Defense_Thoeris type IV,Defense_Tiamat,Defense_Tmn,Defense_Toutatis,Defense_UG1,Defense_UG13a,Defense_UG13b,Defense_UG14,Defense_UG17,Defense_UG19_large,Defense_UG20,Defense_UG21,Defense_UG22,Defense_UG23,Defense_UG29,Defense_UG31,Defense_UG32,Defense_UG33,Defense_UG36,Defense_UG39,Defense_UG4,Defense_UG5_small,Defense_UG6,Defense_UG7,Defense_UG8,Defense_UG9_large,Defense_UG9_small,Defense_Upx,Defense_Uzume,Defense_VP1796,Defense_VP1817,Defense_VP1826,Defense_VP1839,Defense_VP1840,Defense_VP1848,Defense_VP1851,Defense_VP1853,Defense_VSPR,Defense_Viperin solo,Defense_Viperin-TPR,Defense_Wadjet other,Defense_Wadjet type I,Defense_Wadjet type II,Defense_Wadjet type III,Defense_Zorya other,Defense_Zorya type I,Defense_Zorya type II,Defense_Zorya type III,Defense_dCTPdeaminase,Defense_dGTPase,Defense_dXTPase,Defense_gcu142,Defense_gcu233,Defense_gcu24,Defense_gcu76,Defense_gcuWGS21,Defense_pAgo,Defense_pAgo LongA,Defense_pAgo LongB,Defense_pAgo S1A,Defense_pAgo S1B,Defense_pAgo S2B,Defense_pAgo SPARTA,Defense_pAgo other,Defense_pAgo type I,Defense_pAgo type II,Defense_pAgo type III", ",", fixed = TRUE)[[1]],
    name = "Defense"
  ),
  amr = list(
    allowed_cols = strsplit("AMR_aac(2')-IIb,AMR_aac(2')-Ib,AMR_aac(3),AMR_aac(3)-I,AMR_aac(3)-II,AMR_aac(3)-IIIb,AMR_aac(3)-IId,AMR_aac(3)-IV,AMR_aac(6'),AMR_aac(6')-IIa,AMR_aac(6')-Ib,AMR_aac(6')-Ie/aph(2'')-Ia,AMR_aacA8,AMR_aadA1,AMR_aadA16,AMR_aadA2,AMR_aadA27,AMR_aadA5,AMR_aadD1,AMR_aadK,AMR_aadS,AMR_abc-f,AMR_ampC,AMR_ant(2'')-Ia,AMR_ant(3''),AMR_ant(3'')-IIa,AMR_ant(4')-I,AMR_ant(6),AMR_ant(9),AMR_ant(9)-Ic,AMR_aph(3''),AMR_aph(3'')-Ib,AMR_aph(3'),AMR_aph(3')-I,AMR_aph(3')-II,AMR_aph(3')-IIb,AMR_aph(3')-IIc,AMR_aph(3')-Ia,AMR_aph(3')-VIa,AMR_aph(6),AMR_aph(6)-I,AMR_aph(6)-Id,AMR_arr,AMR_arr-3,AMR_bla,AMR_bla-A,AMR_bla-B1-FLAV,AMR_bla2,AMR_bla2a,AMR_blaACT,AMR_blaACT-154,AMR_blaACT-2,AMR_blaADC-18,AMR_blaAIM,AMR_blaAMZ,AMR_blaAQU,AMR_blaASDC,AMR_blaB1PEDO,AMR_blaBPU,AMR_blaCAR,AMR_blaCAU,AMR_blaCESS,AMR_blaCGA,AMR_blaCHM,AMR_blaELC,AMR_blaHER,AMR_blaI,AMR_blaIII,AMR_blaIND,AMR_blaKHM-HMB,AMR_blaL1,AMR_blaL2,AMR_blaLEN-2,AMR_blaLEN-27,AMR_blaLRA1,AMR_blaMAB,AMR_blaMOX,AMR_blaOCH,AMR_blaOXA,AMR_blaOXA-1127,AMR_blaOXA-114g,AMR_blaOXA-272,AMR_blaOXA-815,AMR_blaOXY-4-1,AMR_blaP,AMR_blaPDC-5,AMR_blaPEN-bcc,AMR_blaPJM,AMR_blaPOM,AMR_blaPOM-1,AMR_blaPRC,AMR_blaR1,AMR_blaRATA,AMR_blaSGM,AMR_blaSPR,AMR_blaSRT-3,AMR_blaZ,AMR_cat86,AMR_catA,AMR_catA1,AMR_catB,AMR_catB1,AMR_catB3,AMR_catB7,AMR_catU,AMR_cepH,AMR_cfr,AMR_cipA,AMR_cml,AMR_cmlR,AMR_cmx,AMR_cphA,AMR_cpt,AMR_dfr,AMR_dfrA27,AMR_dfrA33,AMR_dfrA41,AMR_dfrA44,AMR_dfrB,AMR_dfrG,AMR_ere(D),AMR_erm,AMR_erm(30),AMR_erm(31),AMR_erm(A),AMR_erm(B),AMR_erm(C),AMR_erm(O),AMR_erm(T),AMR_erm(V),AMR_estDL136,AMR_floR,AMR_fomB,AMR_fos,AMR_fosA,AMR_fosB,AMR_fosM,AMR_fosX,AMR_fusB,AMR_fusF,AMR_fusH,AMR_gar,AMR_iri,AMR_lnu(A),AMR_lnu(G),AMR_lnu(I),AMR_lsa,AMR_lsa(C),AMR_lsa(D),AMR_mcr-3,AMR_mecA,AMR_mecR1,AMR_mgt,AMR_mph(C),AMR_mph(E),AMR_mphK,AMR_mphM,AMR_mphN,AMR_msr,AMR_msr(A),AMR_msr(E),AMR_mupA,AMR_oqxA,AMR_oqxB,AMR_oqxB9,AMR_otr(A),AMR_qac,AMR_qacA,AMR_qacC,AMR_qacEdelta1,AMR_qacH,AMR_qnrE,AMR_qnrE1,AMR_qnrE3,AMR_qnrE4,AMR_rox,AMR_sal,AMR_satA,AMR_smr,AMR_sul1,AMR_sul2,AMR_sul4,AMR_taeA,AMR_tet,AMR_tet(30),AMR_tet(33),AMR_tet(39),AMR_tet(41),AMR_tet(42),AMR_tet(64),AMR_tet(G),AMR_tet(K),AMR_tet(L),AMR_tet(M),AMR_tet(V),AMR_tetA(58),AMR_tetB(58),AMR_tmexC3,AMR_tmexD,AMR_tmexD3,AMR_toprJ1,AMR_vanA,AMR_vanA-Pa,AMR_vanA-Sc,AMR_vanG,AMR_vanH,AMR_vanH-Sc,AMR_vanJ,AMR_vanK-Sc,AMR_vanR,AMR_vanR-A,AMR_vanR-F,AMR_vanR-O,AMR_vanR-Sc,AMR_vanS,AMR_vanS-Pt,AMR_vanS-Sc,AMR_vanT,AMR_vanW,AMR_vanW-Pt,AMR_vanX,AMR_vanX-Pt,AMR_vanX-Sc,AMR_vanY,AMR_vanY-N,AMR_vanY-Pt,AMR_vanZ-Pt,AMR_vat,AMR_vat(I),AMR_vga,AMR_vga(A),AMR_vgbC,AMR_vmlR,AMR_vph", ",", fixed = TRUE)[[1]],
    name = "AMR"
  ),
  antidefense = list(
    allowed_cols = strsplit("AntiDS_Anti_CBASS,AntiDS_Anti_CRISPR,AntiDS_Anti_Dnd,AntiDS_Anti_Gabija,AntiDS_Anti_Pycsar,AntiDS_Anti_RM,AntiDS_Anti_Thoeris,AntiDS_NADP,AntiDS_Other", ",", fixed = TRUE)[[1]],
    name = "AntiDefense"
  )
)

cat("=== Step 1 (Chromosome+Plasmid overlay): PCoA Calculation by Host and Contig Type ===\n\n")
cat(sprintf("Input file : %s\n", args$input))
cat(sprintf("Gene type  : %s\n", args$type))
if (args$use_counts) {
  cat("Mode       : Using count numbers (Bray-Curtis distance)\n\n")
} else {
  cat("Mode       : Using binary presence/absence (Jaccard distance)\n\n")
}

# ------------------------- 3. Read and filter data -------------------------- #

df <- read.csv(args$input, stringsAsFactors = FALSE, check.names = FALSE)

input_dir <- dirname(normalizePath(args$input, mustWork = TRUE))
output_prefix_base <- basename(args$output)

coord_dir <- file.path(input_dir, "plas_chro_coordinates")
var_dir <- file.path(input_dir, "plas_chro_variance")
for (dir_path in list(coord_dir, var_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", dir_path))
  }
}

if (!"Contig_Type2" %in% colnames(df)) {
  stop("Column 'Contig_Type2' not found in input file.")
}

df_filtered <- df %>%
  filter(Contig_Type2 %in% c("Plasmid", "Chromosome"))

cat(sprintf("Total rows after filtering Contig_Type2 (Plasmid/Chromosome): %d\n\n", nrow(df_filtered)))

# ------------------------- 4. Helper: build matrix -------------------------- #

create_sample_contig_matrix <- function(data, gene_cols) {
  agg_data <- data %>%
    group_by(Sample_ID, Contig_Type2) %>%
    summarise(
      across(all_of(gene_cols), \(x) sum(x, na.rm = TRUE)),
      Host = first(Host),
      Class_CRBC = first(Class_CRBC),
      Order_CRBC = first(Order_CRBC),
      Family_CRBC = first(Family_CRBC),
      Genus_CRBC = first(Genus_CRBC),
      .groups = "drop"
    )

  if (nrow(agg_data) == 0) {
    return(NULL)
  }

  gene_mat <- as.matrix(agg_data[, gene_cols, drop = FALSE])
  rownames(gene_mat) <- paste(agg_data$Sample_ID, agg_data$Contig_Type2, sep = "|")

  meta <- agg_data %>%
    transmute(
      Row_ID = paste(Sample_ID, Contig_Type2, sep = "|"),
      Sample_ID = Sample_ID,
      Contig_Type2 = Contig_Type2,
      Host = Host,
      Class_CRBC = Class_CRBC,
      Order_CRBC = Order_CRBC,
      Family_CRBC = Family_CRBC,
      Genus_CRBC = Genus_CRBC
    )

  list(matrix = gene_mat, metadata = meta)
}

# ------------------------- 5. Single gene type ------------------------------ #

all_cols <- colnames(df_filtered)
cfg <- type_config[[args$type]]
type_name <- args$type

cat("------------------------------------------------------------\n")
cat(sprintf("Processing gene type: %s (%s)\n", cfg$name, type_name))
cat("------------------------------------------------------------\n")

gene_allowed <- intersect(cfg$allowed_cols, all_cols)
missing_allowed <- setdiff(cfg$allowed_cols, all_cols)
if (length(gene_allowed) == 0) {
  stop("No allowed columns for the selected gene type were found in the input file.")
}
if (length(missing_allowed) > 0) {
  cat(sprintf("Warning: %d allowed columns are missing in input and will be skipped.\n", length(missing_allowed)))
}

n_check <- min(100, nrow(df_filtered))
numeric_flags <- sapply(df_filtered[1:n_check, gene_allowed, drop = FALSE], is.numeric)
gene_filtered <- gene_allowed[numeric_flags]

if (length(gene_filtered) == 0) {
  stop("No numeric gene columns were found for the selected gene type.")
}

cat(sprintf("Gene columns (numeric, after filtering): %d\n", length(gene_filtered)))
if (length(gene_filtered) < 2) {
  stop("Fewer than 2 usable gene columns for the selected gene type.")
}

mat_res <- create_sample_contig_matrix(df_filtered, gene_filtered)
if (is.null(mat_res)) {
  stop("No data after aggregation by Sample_ID and Contig_Type2.")
}

gene_mat <- mat_res$matrix
metadata <- mat_res$metadata

row_sums <- rowSums(gene_mat)
gene_mat <- gene_mat[row_sums > 0, , drop = FALSE]
metadata <- metadata[metadata$Row_ID %in% rownames(gene_mat), , drop = FALSE]
cat(sprintf("Rows with at least one gene: %d\n", nrow(gene_mat)))
if (nrow(gene_mat) < 3) {
  stop("Not enough rows for PCoA (need >= 3).")
}

col_sums <- colSums(gene_mat)
gene_mat <- gene_mat[, col_sums > 0, drop = FALSE]
cat(sprintf("Gene types present (non-zero): %d\n", ncol(gene_mat)))
if (ncol(gene_mat) < 2) {
  stop("Not enough gene types for PCoA (need >= 2).")
}

if (args$use_counts) {
  cat("Using count numbers for distance (Bray-Curtis)...\n")
  dist_matrix <- vegdist(gene_mat, method = "bray")
} else {
  cat("Using presence/absence for distance (Jaccard)...\n")
  mat_binary <- (gene_mat > 0) * 1
  dist_matrix <- vegdist(mat_binary, method = "jaccard", binary = TRUE)
}

cat("Running PCoA (Cailliez correction)...\n")
set.seed(123)
pcoa_res <- tryCatch(
  pcoa(dist_matrix, correction = "cailliez"),
  error = function(e) stop(sprintf("Error in PCoA for %s: %s", cfg$name, e$message))
)

eigenvalues <- pcoa_res$values$Eigenvalues
eigenvalues[eigenvalues < 0] <- 0
var_explained <- eigenvalues / sum(eigenvalues) * 100

if (length(var_explained) >= 2) {
  cat(sprintf("Axis1 variance explained: %.2f%%\n", var_explained[1]))
  cat(sprintf("Axis2 variance explained: %.2f%%\n", var_explained[2]))
} else {
  cat("Only one axis available from PCoA; results may be limited.\n")
}

n_axes <- min(2, ncol(pcoa_res$vectors))
if (n_axes < 2) {
  stop("Less than 2 axes available; cannot produce 2D PCoA coordinates.")
}

coords <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
colnames(coords) <- c("PCoA1", "PCoA2")
coords$Row_ID <- rownames(coords)
plot_data <- metadata %>% right_join(coords, by = "Row_ID")

coord_file <- file.path(coord_dir, sprintf("%s_%s_coordinates.csv", output_prefix_base, type_name))
write.csv(plot_data, coord_file, row.names = FALSE)
cat(sprintf("Coordinates saved to: %s\n", coord_file))

var_df <- data.frame(
  Axis = paste0("PCoA", seq_len(min(10, length(var_explained)))),
  Eigenvalue = eigenvalues[seq_len(min(10, length(eigenvalues)))],
  Variance_Explained_Pct = var_explained[seq_len(min(10, length(var_explained)))],
  Cumulative_Variance_Pct = cumsum(var_explained)[seq_len(min(10, length(var_explained)))]
)
var_file <- file.path(var_dir, sprintf("%s_%s_variance.csv", output_prefix_base, type_name))
write.csv(var_df, var_file, row.names = FALSE)
cat(sprintf("Variance explained saved to: %s\n", var_file))
cat("\n=================================================================\n")
cat("PCoA Calculation (chromosome+plasmid overlay design) Complete for selected gene type.\n")
cat("=================================================================\n\n")
cat(sprintf(
  "Gene type: %s  |  Rows: %d  |  Axis1: %.2f%%  |  Axis2: %.2f%%\n",
  type_name, nrow(plot_data), var_explained[1],
  ifelse(length(var_explained) >= 2, var_explained[2], NA_real_)
))
cat("\nOutput directories:\n")
cat(sprintf("  Coordinates: %s\n", coord_dir))
cat(sprintf("  Variance   : %s\n", var_dir))
cat("\nYou can now run '214_step2_pcoa_chro_plas.R' for this gene type (and repeat for the others).\n")