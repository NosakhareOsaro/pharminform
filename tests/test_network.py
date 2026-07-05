from pharminform.network import build_drug_gene_graph, multi_gene_drugs, plot_drug_gene_graph


def test_build_graph_has_gene_and_drug_nodes():
    graph = build_drug_gene_graph()
    kinds = {n: a["kind"] for n, a in graph.nodes(data=True)}
    assert kinds["CYP2D6"] == "gene"
    assert kinds["amitriptyline"] == "drug"


def test_graph_edge_carries_weight_and_directionality():
    graph = build_drug_gene_graph()
    edge = graph["CYP2D6"]["codeine"]
    assert edge["weight"] == 0.10
    assert edge["directionality"] == "prodrug_activation"


def test_multi_gene_drugs_finds_amitriptyline():
    graph = build_drug_gene_graph()
    assert "amitriptyline" in multi_gene_drugs(graph)
    assert "codeine" not in multi_gene_drugs(graph)


def test_plot_drug_gene_graph_writes_file(tmp_path):
    graph = build_drug_gene_graph()
    out = tmp_path / "network.png"
    result_path = plot_drug_gene_graph(graph, out)
    assert result_path.exists()
    assert result_path.stat().st_size > 0
