#!/usr/bin/env Rscript
# 1) Plot isolate network with FR layout
# 2) Styles: compact_disk and scattered
# 3) Optional: center connected nodes when isolated nodes are included
# 4) Optional: resolve node overlaps after layout
# 5) NEW: Shape by Host column in nodes file
# 6) NEW: Annotate nodes with nodes$label text

suppressPackageStartupMessages({
  library(igraph)
  library(qgraph)
  library(RColorBrewer)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_network_isolate_igraph_v2_hostshape_label.R -e edges.tsv -n nodes.tsv -o output.png [options]")
}

# -----------------------------
# Defaults
# -----------------------------
edges_file <- NULL
nodes_file <- NULL
output_file <- NULL

width <- 4000
height <- 4000

edge_width <-1.5
edge_alpha <- 0.8
edge_curved <- 0.2

layout_method <- "igraph"  # qgraph or igraph
niter <- 2000
start_temp_factor <- 2.5

show_labels <- FALSE
title_text <- NULL

color_by <- "origin_class" # origin_class or community
include_isolated <- FALSE

seed <- 1
style <- "compact_disk"    # compact_disk or scattered
vertex_size_by <- "degree" # fixed or degree or strength
vertex_size_fixed <- 2
vertex_size_min <- 2
vertex_size_max <- 5

layout_area_coef <- 20
layout_area_exp <- 1.0
layout_repulse_exp <- 3.2
avoid_overlap <- TRUE
overlap_padding <- 0.1

center_connected <- TRUE
connected_radius <- 0.90
isolated_radius_min <- 0.88
isolated_radius_max <- 0.98

# -----------------------------
# NEW: shape and annotation controls
# -----------------------------
shape_by_host <- TRUE
host_col <- "Host"
host_unknown <- "Unknown"
host_legend_title <- "Host"

annotate_by_label <- TRUE
label_col <- "label"
label_cex <- 0.45
label_color <- "black"
label_font <- 0.2
label_offset_factor <- 0.012

# Base R pch pool
pch_pool <- c(16, 17, 15, 18, 3, 4, 8, 0, 1, 2, 5, 6, 7, 9, 10, 11, 12, 13, 14)

# -----------------------------
# Argument parsing
# -----------------------------
i <- 1
while (i <= length(args)) {
  if (args[i] == "-e" || args[i] == "--edges") {
    edges_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "-n" || args[i] == "--nodes") {
    nodes_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "-o" || args[i] == "--output") {
    output_file <- args[i + 1]; i <- i + 2

  } else if (args[i] == "--width") {
    width <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--height") {
    height <- as.integer(args[i + 1]); i <- i + 2

  } else if (args[i] == "--edge-width") {
    edge_width <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--edge-alpha") {
    edge_alpha <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--edge-curved") {
    edge_curved <- as.numeric(args[i + 1]); i <- i + 2

  } else if (args[i] == "--layout") {
    layout_method <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--niter") {
    niter <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--start-temp-factor") {
    start_temp_factor <- as.numeric(args[i + 1]); i <- i + 2

  } else if (args[i] == "--show-labels") {
    show_labels <- TRUE; i <- i + 1
  } else if (args[i] == "--title") {
    title_text <- args[i + 1]; i <- i + 2

  } else if (args[i] == "--color-by") {
    color_by <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--include-isolated") {
    include_isolated <- TRUE; i <- i + 1
  } else if (args[i] == "--exclude-isolated") {
    include_isolated <- FALSE; i <- i + 1

  } else if (args[i] == "--layout-area-coef") {
    layout_area_coef <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--layout-area-exp") {
    layout_area_exp <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--layout-repulse-exp") {
    layout_repulse_exp <- as.numeric(args[i + 1]); i <- i + 2

  } else if (args[i] == "--seed") {
    seed <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--style") {
    style <- args[i + 1]; i <- i + 2

  } else if (args[i] == "--vertex-size-by") {
    vertex_size_by <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--vertex-size") {
    vertex_size_fixed <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--vertex-size-min") {
    vertex_size_min <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--vertex-size-max") {
    vertex_size_max <- as.numeric(args[i + 1]); i <- i + 2

  } else if (args[i] == "--avoid-overlap") {
    avoid_overlap <- TRUE; i <- i + 1
  } else if (args[i] == "--no-avoid-overlap") {
    avoid_overlap <- FALSE; i <- i + 1
  } else if (args[i] == "--overlap-padding") {
    overlap_padding <- as.numeric(args[i + 1]); i <- i + 2

  } else if (args[i] == "--center-connected") {
    center_connected <- TRUE; i <- i + 1
  } else if (args[i] == "--no-center-connected") {
    center_connected <- FALSE; i <- i + 1
  } else if (args[i] == "--connected-radius") {
    connected_radius <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--isolated-radius-min") {
    isolated_radius_min <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--isolated-radius-max") {
    isolated_radius_max <- as.numeric(args[i + 1]); i <- i + 2

  # NEW: Host shapes
  } else if (args[i] == "--shape-by-host") {
    shape_by_host <- TRUE; i <- i + 1
  } else if (args[i] == "--no-shape-by-host") {
    shape_by_host <- FALSE; i <- i + 1
  } else if (args[i] == "--host-col") {
    host_col <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--host-unknown") {
    host_unknown <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--host-legend-title") {
    host_legend_title <- args[i + 1]; i <- i + 2

  # NEW: annotate labels
  } else if (args[i] == "--annotate-label") {
    annotate_by_label <- TRUE; i <- i + 1
  } else if (args[i] == "--no-annotate-label") {
    annotate_by_label <- FALSE; i <- i + 1
  } else if (args[i] == "--label-col") {
    label_col <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-cex") {
    label_cex <- as.numeric(args[i + 1]); i <- i + 2
  } else if (args[i] == "--label-color") {
    label_color <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--label-font") {
    label_font <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--label-offset-factor") {
    label_offset_factor <- as.numeric(args[i + 1]); i <- i + 2

  } else {
    i <- i + 1
  }
}

if (is.null(edges_file) || is.null(nodes_file) || is.null(output_file)) {
  stop("Missing required arguments: -e/--edges, -n/--nodes, -o/--output")
}

if (!(layout_method %in% c("igraph", "qgraph"))) stop("layout must be 'igraph' or 'qgraph'")
if (!(style %in% c("compact_disk", "scattered"))) stop("style must be 'compact_disk' or 'scattered'")
if (!(vertex_size_by %in% c("fixed", "degree", "strength"))) stop("vertex-size-by must be 'fixed' or 'degree' or 'strength'")
if (!(color_by %in% c("origin_class", "community"))) stop("color-by must be 'origin_class' or 'community'")

# -----------------------------
# Helper functions
# -----------------------------
compute_fr_layout <- function(graph, method, niter, start_temp_factor, area_coef, area_exp, repulse_exp) {
  if (vcount(graph) == 0) return(matrix(numeric(0), ncol = 2))
  if (vcount(graph) == 1) return(matrix(c(0, 0), nrow = 1))
  if (method == "qgraph") {
    e_mat <- as_edgelist(graph, names = FALSE)
    qgraph.layout.fruchtermanreingold(
      e_mat,
      vcount = vcount(graph),
      area = area_coef * (vcount(graph)^area_exp),
      repulse.rad = (vcount(graph)^repulse_exp)
    )
  } else {
    layout_with_fr(
      graph,
      niter = niter,
      start.temp = sqrt(vcount(graph)) * start_temp_factor
    )
  }
}

normalize_layout <- function(layout_xy, mode = "compact_disk", scale_scattered = 1.15) {
  if (nrow(layout_xy) == 0) return(layout_xy)
  layout_xy <- norm_coords(layout_xy, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  if (mode == "compact_disk") {
    r <- sqrt(rowSums(layout_xy^2))
    if (max(r) > 0) layout_xy <- layout_xy / max(r)
    layout_xy <- layout_xy * 0.95
  } else {
    layout_xy <- layout_xy * scale_scattered
  }
  layout_xy
}

resolve_overlaps <- function(layout_xy, radii, padding = 0.4, n_iter = 200, step = 0.02) {
  if (nrow(layout_xy) <= 1) return(layout_xy)
  radii2 <- radii + padding
  n <- nrow(layout_xy)
  for (it in seq_len(n_iter)) {
    moved <- FALSE
    for (ii in 1:(n - 1)) {
      xi <- layout_xy[ii, 1]; yi <- layout_xy[ii, 2]
      for (jj in (ii + 1):n) {
        dx <- layout_xy[jj, 1] - xi
        dy <- layout_xy[jj, 2] - yi
        dist <- sqrt(dx * dx + dy * dy) + 1e-12
        min_dist <- radii2[ii] + radii2[jj]
        if (dist < min_dist) {
          overlap <- (min_dist - dist)
          ux <- dx / dist
          uy <- dy / dist
          push <- overlap * 0.5 * step
          layout_xy[ii, 1] <- layout_xy[ii, 1] - ux * push
          layout_xy[ii, 2] <- layout_xy[ii, 2] - uy * push
          layout_xy[jj, 1] <- layout_xy[jj, 1] + ux * push
          layout_xy[jj, 2] <- layout_xy[jj, 2] + uy * push
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  layout_xy
}

place_isolated_on_ring <- function(n_iso, r_min = 0.88, r_max = 0.98) {
  if (n_iso <= 0) return(matrix(numeric(0), ncol = 2))
  if (n_iso == 1) return(matrix(c(r_max, 0), nrow = 1))
  angles <- seq(0, 2 * pi, length.out = n_iso + 1)[1:n_iso]
  radii <- seq(r_min, r_max, length.out = n_iso)
  cbind(radii * cos(angles), radii * sin(angles))
}

# -----------------------------
# Step 1: Load data
# -----------------------------
cat("Loading edges and nodes data...\n")
edges_df <- read.delim(edges_file, sep = "\t", stringsAsFactors = FALSE)
nodes_df <- read.delim(nodes_file, sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d edges and %d nodes\n", nrow(edges_df), nrow(nodes_df)))

# -----------------------------
# Step 2: Normalize edge columns and build graph
# -----------------------------
cat("Creating igraph graph...\n")
if ("g1" %in% colnames(edges_df) && "g2" %in% colnames(edges_df)) {
  edges_df$source <- edges_df$g1
  edges_df$target <- edges_df$g2
} else if ("source" %in% colnames(edges_df) && "target" %in% colnames(edges_df)) {
  # keep
} else {
  stop("Edge file must have g1/g2 or source/target columns")
}

if (!("id" %in% colnames(nodes_df))) stop("Nodes file must have 'id' column")

edge_node_ids <- unique(c(edges_df$source, edges_df$target))
if (include_isolated) {
  cat("Including isolated nodes...\n")
  nodes_df_filtered <- nodes_df
} else {
  cat("Excluding isolated nodes...\n")
  nodes_df_filtered <- nodes_df[nodes_df$id %in% edge_node_ids, ]
}
nodes_df_filtered <- nodes_df_filtered[!duplicated(nodes_df_filtered$id), ]

g <- graph_from_data_frame(
  d = edges_df[, c("source", "target")],
  vertices = nodes_df_filtered,
  directed = FALSE
)

if ("jaccard" %in% colnames(edges_df)) {
  E(g)$weight <- edges_df$jaccard
  cat("Using jaccard similarity as edge weights\n")
} else {
  E(g)$weight <- NA_real_
  cat("No edge weights found, using unweighted graph\n")
}

deg_all <- degree(g)
isolated_count <- sum(deg_all == 0)
cat(sprintf("Created graph with %d vertices and %d edges\n", vcount(g), ecount(g)))
if (isolated_count > 0) cat(sprintf("Found %d isolated nodes\n", isolated_count))

# -----------------------------
# Step 3: Colors
# -----------------------------
vertex_colors <- rep("#808080", vcount(g))
edge.col <- rep("#999999", ecount(g))

if (color_by == "origin_class") {
  cat("Using origin_class for node colors...\n")

  class_colors <- c(
    "Actinomycetia" = "#98df8a",
    "Alphaproteobacteria" = "#aec7e8",
    "Bacilli" = "#ff7f0e",
    "Thermoleophilia" = "#ff9896",
    "Bacteroidia" = "#d62728",
    "Gammaproteobacteria" = "#ffbb78",
    "Deinococci" = "#1f77b4",
    "Acidimicrobiia" = "#2ca02c",
    "Campylobacteria" = "#9467bd",
    "Unknown" = "#808080"
  )

  if ("origin_class" %in% colnames(nodes_df_filtered)) {
    node_cls <- nodes_df_filtered$origin_class[match(V(g)$name, nodes_df_filtered$id)]
  } else if ("Class_CRBC" %in% colnames(nodes_df_filtered)) {
    node_cls <- nodes_df_filtered$Class_CRBC[match(V(g)$name, nodes_df_filtered$id)]
  } else if ("Class" %in% colnames(nodes_df_filtered)) {
    node_cls <- nodes_df_filtered$Class[match(V(g)$name, nodes_df_filtered$id)]
  } else {
    node_cls <- rep("Unknown", vcount(g))
  }
  node_cls[is.na(node_cls) | node_cls == ""] <- "Unknown"

  vertex_colors <- ifelse(node_cls %in% names(class_colors), class_colors[node_cls], "#808080")

  ee <- ends(g, es = E(g), names = FALSE)
  edge.col <- vertex_colors[ee[, 1]]

} else {
  cat("Detecting communities using Louvain...\n")
  comm <- cluster_louvain(g)
  mem <- membership(comm)
  k <- max(mem)

  pal <- brewer.pal(min(12, max(3, k)), "Set3")
  if (k > length(pal)) pal <- colorRampPalette(pal)(k)
  vertex_colors <- pal[mem]

  ee <- ends(g, es = E(g), names = FALSE)
  edge.col <- vertex_colors[ee[, 1]]
}

# -----------------------------
# Step 4: Vertex sizes
# -----------------------------
cat(sprintf("Vertex size mode: %s\n", vertex_size_by))

if (vertex_size_by == "fixed") {
  vertex_sizes <- rep(vertex_size_fixed, vcount(g))
} else if (vertex_size_by == "degree") {
  if (max(deg_all) == 0) vertex_sizes <- rep(vertex_size_fixed, vcount(g))
  else vertex_sizes <- rescale(deg_all, to = c(vertex_size_min, vertex_size_max))
} else {
  if (all(is.na(E(g)$weight))) {
    if (max(deg_all) == 0) vertex_sizes <- rep(vertex_size_fixed, vcount(g))
    else vertex_sizes <- rescale(deg_all, to = c(vertex_size_min, vertex_size_max))
  } else {
    st <- strength(g, weights = E(g)$weight)
    if (max(st) == 0) vertex_sizes <- rep(vertex_size_fixed, vcount(g))
    else vertex_sizes <- rescale(st, to = c(vertex_size_min, vertex_size_max))
  }
}

radii <- vertex_sizes / max(vertex_sizes) * 0.06

# -----------------------------
# Step 4.5: Shapes by Host
# -----------------------------
vertex_pch <- rep(16, vcount(g))
host_levels <- character(0)
host_to_pch <- NULL

if (shape_by_host) {
  if (!(host_col %in% colnames(nodes_df_filtered))) {
    stop(sprintf("shape-by-host enabled but nodes file missing column: %s", host_col))
  }
  host_vec <- nodes_df_filtered[[host_col]][match(V(g)$name, nodes_df_filtered$id)]
  host_vec[is.na(host_vec) | host_vec == ""] <- host_unknown

  host_levels <- sort(unique(host_vec))
  host_to_pch <- setNames(pch_pool[((seq_along(host_levels) - 1) %% length(pch_pool)) + 1], host_levels)
  vertex_pch <- unname(host_to_pch[host_vec])
}

# -----------------------------
# Step 4.6: Annotation by label
# -----------------------------
vertex_annot <- rep(NA_character_, vcount(g))
if (annotate_by_label) {
  if (!(label_col %in% colnames(nodes_df_filtered))) {
    stop(sprintf("annotate-label enabled but nodes file missing column: %s", label_col))
  }
  lab_vec <- nodes_df_filtered[[label_col]][match(V(g)$name, nodes_df_filtered$id)]
  lab_vec[is.na(lab_vec) | lab_vec == ""] <- NA_character_
  vertex_annot <- lab_vec
}

# -----------------------------
# Step 5: Layout
# -----------------------------
cat(sprintf("Setting seed: %d\n", seed))
set.seed(seed)

layout_fr <- NULL

if (include_isolated && center_connected && style == "compact_disk") {
  cat("Center-connected mode enabled for include_isolated + compact_disk\n")

  idx_conn <- which(deg_all > 0)
  idx_iso <- which(deg_all == 0)

  layout_fr <- matrix(0, nrow = vcount(g), ncol = 2)

  if (length(idx_conn) > 0) {
    g_conn <- induced_subgraph(g, idx_conn)
    lay_conn <- compute_fr_layout(g_conn, layout_method, niter, start_temp_factor,
                                 layout_area_coef, layout_area_exp, layout_repulse_exp)
    lay_conn <- normalize_layout(lay_conn, mode = "compact_disk")
    lay_conn <- lay_conn * connected_radius
    layout_fr[idx_conn, ] <- lay_conn
  }

  if (length(idx_iso) > 0) {
    lay_iso <- place_isolated_on_ring(length(idx_iso), isolated_radius_min, isolated_radius_max)
    layout_fr[idx_iso, ] <- lay_iso
  }

} else {
  cat(sprintf("Calculating layout: %s, method: %s\n", style, layout_method))
  layout_fr <- compute_fr_layout(g, layout_method, niter, start_temp_factor,
                                 layout_area_coef, layout_area_exp, layout_repulse_exp)
  layout_fr <- normalize_layout(layout_fr, mode = style)
}

if (avoid_overlap && style == "compact_disk") {
  cat("Resolving overlaps...\n")
  layout_fr <- resolve_overlaps(layout_fr, radii = radii, padding = overlap_padding, n_iter = 180, step = 0.03)
  layout_fr <- normalize_layout(layout_fr, mode = "compact_disk")
}

# -----------------------------
# Step 6: Output
# -----------------------------
cat(sprintf("Plotting network to %s...\n", output_file))

output_format <- tools::file_ext(output_file)
if (output_format == "png") {
  png(output_file, width = width, height = height, res = 300)
} else if (output_format == "pdf") {
  pdf(output_file, width = width / 100, height = height / 100)
} else if (output_format == "svg") {
  svg(output_file, width = width / 100, height = height / 100)
} else {
  png(output_file, width = width, height = height, res = 300)
}

plot_title <- if (is.null(title_text)) "Isolate Network Graph" else title_text

par(bg = "white", mar = c(2, 2, 3, 2))

if (style == "compact_disk") {
  xlim_use <- c(-1, 1)
  ylim_use <- c(-1, 1)
  rescale_use <- FALSE
  asp_use <- 1
} else {
  xlim_use <- NULL
  ylim_use <- NULL
  rescale_use <- TRUE
  asp_use <- 0
}

# -----------------------------
# Plot edges first, then points and text
# Key: vertex.size = 0 prevents edge endpoint trimming
# -----------------------------
plot.igraph(
  g,
  layout = layout_fr,
  rescale = rescale_use,
  xlim = xlim_use,
  ylim = ylim_use,
  asp = asp_use,
  vertex.shape = "none",
  vertex.size = 0,
  vertex.label = NA,
  edge.color = adjustcolor(edge.col, alpha.f = edge_alpha),
  edge.width = edge_width,
  edge.curved = edge_curved,
  main = plot_title,
  cex.main = 1.0,
  font.main = 2
)

# Draw nodes with Host-dependent pch
points(
  layout_fr[, 1], layout_fr[, 2],
  pch = vertex_pch,
  cex = vertex_sizes / 2.5,
  col = vertex_colors
)

# Optional label annotation
if (annotate_by_label && any(!is.na(vertex_annot))) {
  dx <- diff(par("usr")[1:2])
  dy <- diff(par("usr")[3:4])
  text(
    layout_fr[, 1],
    layout_fr[, 2] + (label_offset_factor * max(dx, dy)),
    labels = vertex_annot,
    cex = label_cex,
    col = label_color,
    font = label_font
  )
}

# Host legend (shapes)
if (shape_by_host && length(host_levels) > 0) {
  legend(
    "topleft",
    legend = host_levels,
    pch = unname(host_to_pch[host_levels]),
    col = "black",
    pt.cex = 1.1,
    cex = 0.7,
    bty = "n",
    title = host_legend_title,
    inset = c(0.02, 0.02)
  )
}

dev.off()

cat("Network visualization completed successfully\n")
cat(sprintf("Output saved to: %s\n", output_file))
cat(sprintf("Host shapes: %s (column: %s)\n", ifelse(shape_by_host, "ON", "OFF"), host_col))
cat(sprintf("Label annotation: %s (column: %s)\n", ifelse(annotate_by_label, "ON", "OFF"), label_col))
