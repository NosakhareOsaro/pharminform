"""NetworkX drug-gene interaction network.

Builds a bipartite graph (genes <-> drugs) from
data/drug_pathway_weights.json, where an edge's weight is the published
fraction_metabolized (or None if not published, i.e. equal-weight fallback
territory -- see composite_risk.py) and its directionality attribute flags
prodrug-activation pathways where reduced gene function does NOT increase
toxicity risk (see DATA_SOURCES.md).
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx

from .composite_risk import load_drug_pathway_weights


def build_drug_gene_graph(weights_data: dict | None = None) -> nx.DiGraph:
    data = weights_data if weights_data is not None else load_drug_pathway_weights()
    graph = nx.DiGraph()
    for entry in data["drugs"]:
        drug = entry["drug"]
        graph.add_node(drug, kind="drug", rxnorm_id=entry.get("rxnorm_id"))
        for gene_info in entry["contributing_genes"]:
            gene = gene_info["gene"]
            graph.add_node(gene, kind="gene")
            graph.add_edge(
                gene,
                drug,
                weight=gene_info["fraction_metabolized"],
                directionality=entry.get("directionality", "clearance"),
            )
    return graph


def multi_gene_drugs(graph: nx.DiGraph) -> list[str]:
    """Drugs with more than one contributing gene -- the composite model's
    core use case."""
    return [
        node for node, attrs in graph.nodes(data=True)
        if attrs.get("kind") == "drug" and graph.in_degree(node) > 1
    ]


def plot_drug_gene_graph(graph: nx.DiGraph, out_path: str | Path) -> Path:
    """Render the graph to a static image file (PNG/SVG/PDF, by extension).
    Requires matplotlib; genes and drugs are colored differently, and
    prodrug-activation edges are drawn dashed to visually distinguish them
    from clearance-pathway edges.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gene_nodes = [n for n, a in graph.nodes(data=True) if a.get("kind") == "gene"]
    drug_nodes = [n for n, a in graph.nodes(data=True) if a.get("kind") == "drug"]

    pos = nx.bipartite_layout(graph, gene_nodes) if gene_nodes and drug_nodes else nx.spring_layout(graph)

    fig, ax = plt.subplots(figsize=(10, 8))
    nx.draw_networkx_nodes(graph, pos, nodelist=gene_nodes, node_color="#4C72B0", node_shape="s", node_size=1200, ax=ax)
    nx.draw_networkx_nodes(graph, pos, nodelist=drug_nodes, node_color="#DD8452", node_size=900, ax=ax)

    solid_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get("directionality") != "prodrug_activation"]
    dashed_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get("directionality") == "prodrug_activation"]
    nx.draw_networkx_edges(graph, pos, edgelist=solid_edges, ax=ax, arrows=True)
    nx.draw_networkx_edges(graph, pos, edgelist=dashed_edges, style="dashed", ax=ax, arrows=True)

    nx.draw_networkx_labels(graph, pos, font_size=8, ax=ax)
    ax.set_title("PharmInform: gene → drug pathway network\n(dashed = prodrug activation, see DATA_SOURCES.md)")
    ax.axis("off")

    out_path = Path(out_path)
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return out_path
