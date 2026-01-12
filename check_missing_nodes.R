#!/usr/bin/env Rscript
edges <- read.delim("Result/NCBI_4395_Batch/07_Network/isolate_network/Annotation/All_genome_sharing_edges.tsv", sep="\t", stringsAsFactors=FALSE)
nodes <- read.delim("Result/NCBI_4395_Batch/07_Network/isolate_network/Annotation/All_genome_sharing_nodes.tsv", sep="\t", stringsAsFactors=FALSE)

edge_nodes <- unique(c(edges$g1, edges$g2))
node_ids <- unique(nodes$id)
missing <- setdiff(edge_nodes, node_ids)

cat(sprintf("Nodes in edges file: %d\n", length(edge_nodes)))
cat(sprintf("Nodes in nodes file: %d\n", length(node_ids)))
cat(sprintf("Nodes in edges but not in nodes: %d\n", length(missing)))

if(length(missing) > 0) {
  cat("Missing node IDs:\n")
  for(i in 1:min(20, length(missing))) {
    cat(sprintf("  %s\n", missing[i]))
  }
  if(length(missing) > 20) {
    cat(sprintf("  ... and %d more\n", length(missing) - 20))
  }
}
