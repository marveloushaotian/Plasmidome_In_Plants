#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from itertools import combinations

import pandas as pd

try:
    import networkx as nx
except ImportError:
    nx = None

GENOME_RE = re.compile(r"^(CRBC_G\d+)")


def extract_genome_id(contig_id: str) -> str:
    s = str(contig_id).strip()
    m = GENOME_RE.match(s)
    if m:
        return m.group(1)
    # fallback: split on first '_' if needed
    return s.split("_")[0]


def load_cluster_map(path: str) -> pd.DataFrame:
    """
    Reads 2-column mapping:
      col0 = cluster representative (cluster ID)
      col1 = member contig ID
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] < 2:
        raise ValueError(f"{path} must have >=2 columns.")
    df = df.iloc[:, :2].copy()
    df.columns = ["cluster", "contig_id"]
    df["cluster"] = df["cluster"].astype(str).str.strip()
    df["contig_id"] = df["contig_id"].astype(str).str.strip()
    return df


def load_meta(meta_path: str) -> pd.DataFrame:
    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    cols = {c.lower(): c for c in meta.columns}
    if "genomeid_standard" in cols:
        id_col = cols["genomeid_standard"]
    else:
        id_col = meta.columns[0]
    meta = meta.set_index(id_col)
    meta.index = meta.index.astype(str).str.strip()
    # strip all cells
    for c in meta.columns:
        meta[c] = meta[c].astype(str).str.strip()
    return meta


def build_nodes(cluster_df: pd.DataFrame, meta: pd.DataFrame, tax_cols):
    nodes = pd.DataFrame({"contig": pd.unique(cluster_df["contig_id"])})
    nodes["genome"] = nodes["contig"].apply(extract_genome_id)

    # attach taxonomy
    for col in tax_cols:
        if col in meta.columns:
            nodes[col] = nodes["genome"].map(meta[col]).fillna("NA")
        else:
            nodes[col] = "NA"

    # also record host if present
    if "Host" in meta.columns and "Host" not in tax_cols:
        nodes["Host"] = nodes["genome"].map(meta["Host"]).fillna("NA")

    return nodes


def build_edges(cluster_df: pd.DataFrame, mode="star",
                min_cluster_size=2, max_cluster_size=0):
    """
    mode:
      - star: edges (cluster_rep -> each member) excluding self
      - clique: all pairwise edges among members (undirected) excluding self
    min_cluster_size: ignore clusters smaller than this
    max_cluster_size: if >0, ignore clusters larger than this (safety)
    """
    edges = []

    # group members by cluster
    grouped = cluster_df.groupby("cluster")["contig_id"].apply(list)

    for cl, members in grouped.items():
        # deduplicate while preserving order
        seen = set()
        mem = []
        for x in members:
            if x not in seen:
                mem.append(x)
                seen.add(x)

        n = len(mem)
        if n < min_cluster_size:
            continue
        if max_cluster_size and n > max_cluster_size:
            continue

        if mode == "star":
            rep = cl
            for m in mem:
                if m == rep:
                    continue
                edges.append({"source": rep, "target": m, "cluster": cl})
        else:  # clique
            for a, b in combinations(mem, 2):
                if a == b:
                    continue
                edges.append({"source": a, "target": b, "cluster": cl})

    edf = pd.DataFrame(edges)
    if edf.empty:
        edf = pd.DataFrame(columns=["source", "target", "cluster"])
    return edf


def write_graphml(nodes_df, edges_df, out_graphml, node_cols):
    if nx is None:
        print("⚠️ networkx not installed; skipping GraphML.")
        return

    G = nx.Graph()

    # add nodes with attributes
    for _, r in nodes_df.iterrows():
        attrs = {c: r[c] for c in node_cols if c in nodes_df.columns}
        G.add_node(r["contig"], **attrs)

    # add edges (undirected)
    for _, r in edges_df.iterrows():
        u = r["source"]
        v = r["target"]
        if u == v:
            continue
        # multiple clusters between same contigs is unlikely here; still safe:
        if G.has_edge(u, v):
            # store as list of clusters if repeated
            prev = G[u][v].get("clusters", "")
            if prev:
                G[u][v]["clusters"] = prev + "," + str(r["cluster"])
            else:
                G[u][v]["clusters"] = str(r["cluster"])
            G[u][v]["weight"] = float(G[u][v].get("weight", 1.0)) + 1.0
        else:
            G.add_edge(u, v, clusters=str(r["cluster"]), weight=1.0)

    nx.write_graphml(G, out_graphml)
    print(f"✅ GraphML → {out_graphml} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")


def main():
    ap = argparse.ArgumentParser(description="Build contig-level plasmid linkage network from MMseqs cluster mapping.")
    ap.add_argument("--cluster-map", required=True,
                    help="*_cluster_output_cluster.renamed.tsv (2 cols: cluster_rep, member_contig).")
    ap.add_argument("--meta", required=True,
                    help="meta_statistics.txt (GenomeID_standard + taxonomy columns).")
    ap.add_argument("--tax-cols", default="Kingdom,Phylum,Class,Order,Family,Genus",
                    help="Comma-separated taxonomy columns to attach to nodes.")
    ap.add_argument("--mode", choices=["star", "clique"], default="star",
                    help="Edge mode: star=rep→member edges (lighter), clique=all pairs per cluster (heavier).")
    ap.add_argument("--min-cluster-size", type=int, default=2,
                    help="Ignore clusters with < this many members.")
    ap.add_argument("--max-cluster-size", type=int, default=0,
                    help="If >0, ignore clusters with > this many members (safety for clique).")
    ap.add_argument("--out-prefix", required=True,
                    help="Output prefix (folder + basename).")
    ap.add_argument("--graphml", action="store_true",
                    help="Also write GraphML (requires networkx).")
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_prefix)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    tax_cols = [c.strip() for c in args.tax_cols.split(",") if c.strip()]

    cm = load_cluster_map(args.cluster_map)
    meta = load_meta(args.meta)

    # Nodes (all contigs appearing in mapping, both reps and members)
    # Note: reps are also contigs (often self-mapped), but just in case:
    all_contigs = pd.unique(pd.concat([cm["contig_id"], cm["cluster"]], ignore_index=True))
    cm2 = pd.DataFrame({"cluster": cm["cluster"], "contig_id": cm["contig_id"]})
    # Ensure reps included as contigs too (even if no explicit self line)
    reps = pd.DataFrame({"cluster": cm["cluster"].unique(), "contig_id": cm["cluster"].unique()})
    cm2 = pd.concat([cm2, reps], ignore_index=True).drop_duplicates()

    nodes_df = build_nodes(cm2, meta, tax_cols)
    edges_df = build_edges(cm2, mode=args.mode,
                           min_cluster_size=args.min_cluster_size,
                           max_cluster_size=args.max_cluster_size)

    nodes_out = args.out_prefix + "_contig_nodes.tsv"
    edges_out = args.out_prefix + "_contig_edges.tsv"
    nodes_df.to_csv(nodes_out, sep="\t", index=False)
    edges_df.to_csv(edges_out, sep="\t", index=False)

    print(f"✅ Nodes → {nodes_out} (n={nodes_df.shape[0]})")
    print(f"✅ Edges → {edges_out} (n={edges_df.shape[0]})")

    if args.graphml:
        out_graphml = args.out_prefix + "_contig.graphml"
        node_cols = ["genome"] + tax_cols + (["Host"] if "Host" in nodes_df.columns else [])
        write_graphml(nodes_df, edges_df, out_graphml, node_cols=node_cols)


if __name__ == "__main__":
    main()

