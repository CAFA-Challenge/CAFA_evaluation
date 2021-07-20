from math import log
from typing import Iterable
from goatools.obo_parser import GODag
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def get_graph_root(graph):
    graph_root = [node for node, degree in graph.in_degree() if degree == 0]
    return(graph_root[0])


def get_nodes_ia(nodes: Iterable[str], graph: nx.DiGraph) -> float:
    ia_vals = [
        node_obj.get("ia")
        for term, node_obj in graph.nodes.items()
        if term in nodes
    ]
    return sum(ia_vals)


def get_all_upper_edges(goatools_node, dag: GODag):
    """This is a fork of goatools get_all_parent_edges() that takes into
    account both is_a and part_of relationships where get_all_parent_edges()
    does NOT include part_of relationships.

    Return tuples for all parent GO IDs, containing current GO ID and parent GO ID."""

    all_upper_edges = set()
    node_id = goatools_node.item_id


    #print(f"{node_id}\t{dag[node_id]}")

    #for parent in dag[node_id].get_goterms_upper_rels({"is_a", "part_of"}):
    for parent in dag[node_id].get_goterms_upper_rels({"is_a", "part_of"}):
        parent_id = parent.item_id
        all_upper_edges.add((goatools_node.item_id, parent_id))
        all_upper_edges |= get_all_upper_edges(dag[parent_id], dag)

    return all_upper_edges


def get_annotation_dataframe(
        propagation_map_filepath: str,
        dag: GODag,
        annotation_pairs: Iterable,
        pare_df: bool = False,
) -> pd.DataFrame:
    """Creates and returns a pandas DataFrame representing propagated protein
    annotations.

    the resulting DataFrame has this form:
                       GO:0000123  GO:0000124  GO:0000126  GO:0000137  GO:0000139
    T100900000026           0           0           0           0           0
    T100900000115           0           0           0           0           0
    T100900000116           0           0           0           0           0
    T100900000167           0           0           0           0           0
    T100900000453           0           0           0           0           0
    """
    proteins = sorted(list({protein for protein, _ in annotation_pairs}))

    annotation_df = pd.DataFrame(
        data=0, index=proteins, columns=sorted(list(dag.keys()))
    )
    propagation_map = pd.read_pickle(propagation_map_filepath)
    # propagation_map = pd.read_pickle("./cco_propagation_map.pkl")

    # loop the leaf-only groundtruth, updating the dataframe with propagated
    # annotations (toggle zero vals to 1 where appropriate):
    for protein, node in annotation_pairs:
        # TODO: Consider alt-IDs
        try:
            propagation_df_row = propagation_map.loc[node, :]
            propagation_mask = propagation_df_row == 1
            # column names (GO IDs) that are ancestors of our leaf GO ID:
            propagation_mask_cols = propagation_mask.index[propagation_mask].tolist()
            annotation_df.loc[protein, propagation_mask_cols] = 1
        except KeyError:
            pass

    if pare_df:
        # TODO: For the sake of expediency, temporarily pare the terms down to those
        # TODO: that are actually used by our limited protein dataset:
        term_mask = annotation_df.sum(axis="index") > 0
        cols_of_interest = term_mask.index[term_mask].tolist()
        annotation_df = annotation_df.loc[:, cols_of_interest]

    return annotation_df


'''
def get_graph(annotation_pairs: Iterable, dag: GODag) -> nx.DiGraph:
    """Constructs and returns a nx.DiGraph of the propagated ontology terms
    (terms from annotation pairs param).

    """
    graph = nx.DiGraph()

    # Add nodes and edges to the graph:
    # annotation_pairs should be "leaf only" data...
    for _, leaf_node in annotation_pairs:

        try:
            # ancestry is a set of two-element tuples:
            # {('GO:0044464', 'GO:0005575'), ('GO:0044424', 'GO:0044464'),...
            # where the first item is a child and the second is an immediate parent
            ancestry = get_all_upper_edges(dag[leaf_node], dag)

            # this set will include all nodes INCLUDING the "leaf":
            unique_nodes = {term for pair in ancestry for term in pair}

            # Add a value to aid in rendering the graph "root-down",
            # there is probably a better way to do this:
            for node in unique_nodes:
                graph.add_node(node, weight=-dag[node].depth)

            # the node pairs in ancestry are child->parent,
            # while networkx edges should be the inverse,
            # so account for that when creating the edges:
            graph.add_edges_from({(parent, child) for child, parent in ancestry})

        except KeyError:
            # TODO: see note about GO:0098686 elsewhere
            pass

    return graph
'''


def get_graph(annotation_dataframe: pd.DataFrame, dag: GODag) -> nx.DiGraph:
    """Constructs and returns a nx.DiGraph of the propagated ontology terms
    (terms from annotation pairs param).

    """
    graph = nx.DiGraph()

    # Extract the non-zero ontology terms from the DataFrame:
    nonzero_mask = annotation_dataframe != 0
    # get a list of the annotated column names (GO IDs) base on the mask:
    nonzero_mask = nonzero_mask.any(axis="index")
    annotated_go_ids = annotation_dataframe.columns[nonzero_mask].tolist()

    # Add nodes and edges to the graph:
    # annotation_pairs should be "leaf only" data...
    for node in annotated_go_ids:

        try:
            # ancestry is a set of two-element tuples:
            # {('GO:0044464', 'GO:0005575'), ('GO:0044424', 'GO:0044464'),...
            # where the first item is a child and the second is an immediate parent
            ancestry = get_all_upper_edges(dag[node], dag)

            # this set will include all nodes INCLUDING the "leaf":
            unique_nodes = {term for pair in ancestry for term in pair}

            # Add a value to aid in rendering the graph "root-down",
            # there is probably a better way to do this:
            for node in unique_nodes:
                graph.add_node(node, weight=-dag[node].depth)

            # the node pairs in ancestry are child->parent,
            # while networkx edges should be the inverse,
            # so account for that when creating the edges:
            graph.add_edges_from({(parent, child) for child, parent in ancestry})

        except KeyError:
            # TODO: see note about GO:0098686 elsewhere
            pass

    return graph


def get_parent_terms(
    term_id: str, dag: GODag, relationships: Iterable[str] = []
) -> set:
    """ Returns the parent node names for the given node over multiple relationship types and
    not only "is_a". """
    try:
        parent_terms = {parent.item_id for parent in dag[term_id].parents}
        for relation_key in relationships:
            parent_terms.update(
                {
                    upper.item_id
                    for upper in dag[term_id].relationship.get(relation_key, [])
                }
            )
    except KeyError:
        parent_terms = None

    return parent_terms


def compute_nodes_information_content(
        graph: nx.DiGraph,
        annotation_dataframe: pd.DataFrame,
        dag: GODag
) -> None:
    
    graph_root = [node for node, degree in graph.in_degree() if degree == 0]
    
    for node_label in graph.nodes:
        node = graph.nodes[node_label]

        if node_label in graph_root:
            node["precision"] = 0
            node["ia"] = 0
        else:
            parent_terms = get_parent_terms(
                term_id=node_label, dag=dag, relationships=("part_of",)
            )
            # this is a boolean matrix of all proteins x the GO terms of interest
            # (terms that are parents of the term being evaluated:
            #print(parent_terms)
            parent_mask = annotation_dataframe.loc[:, parent_terms].all(axis=1)
            parent_count = len(annotation_dataframe.loc[parent_mask]) + 1
            term_count = sum(annotation_dataframe.loc[:, node_label]) + 1
            node["precision"] = term_count/parent_count
            node["ia"] = -log(node["precision"])


def compute_nodes_information_content_BAK(
        graph: nx.DiGraph, annotation_dataframe: pd.DataFrame
) -> None:
    """ Decorates each node of the graph with data based on annotation_dataframe

    Parameters
    ---------
    graph
        graph representing the ontological hierarchy of terms used for annotation
    annotation_dataframe
        dataframe representing propagated annotation of a set of proteins

    Returns
    ---------
    None
    """
    # Here we compute Pr(n2|n1) and ia (-log2) for each node in the graph:
    for node_label in graph.nodes:
        node = graph.nodes[node_label]
        # Compute the Pr(d|bc) value:

        # count the occurences of the node:
        node_annotation_sum = sum(annotation_dataframe.loc[:, node_label]) + 1

        # What are the direct ancestors for the node at hand?
        node_all_ancestors = [
            parent for parent, child in graph.edges if child == node_label
        ]

        # TODO: Review this comment for accuracy and review the code too.

        # Take a slice of the dataframe that is ALL proteins X ancestor GO IDs
        # associated with the current node of interest. It will have this form:
        #
        #                       GO:0044444  GO:0044430  GO:0005815
        # T100900000026           0           0           0
        # T100900000115           1           0           0
        # T100900000116           0           0           0
        # T100900000167           1           0           0
        # T100900000453           0           0           0
        # ...                   ...         ...         ...
        #
        # We do this to calculate "probability" which is the ratio of the
        # sum of annotations for the term of interest (TOI) over
        # the sum of annotations for the TOI's direct parent terms
        # See figure 1 in the paper
        # @ https://academic.oup.com/bioinformatics/article/29/13/i53/195366
        # for more detail.

        if len(node_all_ancestors) == 0:
            # this is a root node for the namespace:
            precision = 1
        else:
            node_all_ancestors_df = annotation_dataframe.loc[:, node_all_ancestors]

            # Count the rows where ALL direct ancestors == 1:
            ancestor_annotation_sum = sum(node_all_ancestors_df.all(axis="columns"))

            ancestor_annotation_sum += len(node_all_ancestors)

            precision = round(node_annotation_sum / ancestor_annotation_sum, 3)

        node["pr"] = precision
        # Finally, we add a field for -log2 which is used as a matter of convention:
        # node["nl2"] = round(-log2(node["pr"]), 3)


def compute_protein_information_content(
        protein_id: str, graph: nx.DiGraph, annotation_dataframe: pd.DataFrame
) -> float:
    """Aggregate and return the information content for each node associated
    with the given protein
    """
    try:
        # The annotation_dataframe param maps proteins (row indices) to
        # GO IDs (columns). This DataFrame should contain binary data
        # representing annotation (1) vs NO annotation (0).
        # However, this function should work for other values as long as the
        # no annotation state is represented with 0.

        # Create a boolean pandas.Series for a single protein
        # indicating annotated/not-annotated:
        nonzero_mask = annotation_dataframe.loc[protein_id] != 0
        # get a list of the annotated column names (GO IDs) base on the mask:
        annotated_go_ids = annotation_dataframe.columns[nonzero_mask].tolist()

        # finally use the annotated GO IDs to read data (nl2 values) from
        # the relevant propagated graph nodes.
        # information content for a protein is the "accretion" of nl2 values:
        info_content = sum(
            [-log(graph.nodes[node].get("pr")) for node in graph if node in annotated_go_ids if
             graph.nodes[node].get("pr") != 0]
        )

        return info_content

    except KeyError:
        # TODO:
        # For CAFA3 CCO "leaf only" ground truth, there is an odd case of
        # GO:0098686 that is NOT in the corresponding obo file. That
        # case (and probably others) are carelessly squashed here:
        return 0


def render_graph(graph: nx.DiGraph) -> None:
    """Stub function for rendering a networkx DiGraph.
    Needs to be fleshed out.

    Graph nodes are expected to have certain fields:
    1. a 'weight' value to aid in rendering
    2. pr (calculated probability, see compute_nodes_information_content())
    3. nl2 (negative log2, see compute_nodes_information_content())
    """

    node_colors = [graph.nodes[n].get("color", "gray") for n in graph.nodes]

    pos = nx.multipartite_layout(graph, subset_key="weight", align="horizontal")
    nx.draw_networkx_nodes(
        graph, pos, cmap=plt.get_cmap("jet"), node_color=node_colors, node_size=400
    )
    pos2 = {k: [v[0], v[1]] for k, v in pos.items()}
    nx.draw_networkx_labels(graph, pos2, font_size=11, font_color="black")
    nx.draw_networkx_edges(graph, pos, edge_color="gray", arrowsize=20, arrows=True)
    # nx.draw_networkx_labels(graph, pos, font_size=12, font_color="black")
    # nx.draw_networkx_edges(graph, pos, edge_color="gray", arrows=True)

    # Add labels to the graph for the Pr() values:
    precision_labels = {k: graph.nodes[k].get("pr", 0) for k in graph}
    precision_labels_pos = {k: (v[0], v[1] - 0.05) for k, v in pos.items()}

    nx.draw_networkx_labels(
        graph, precision_labels_pos, precision_labels, font_color="maroon", font_size=10
    )
    # Add labels to the graph for the negative log2 values:

    #neg_log2_labels = {k: round(-log2(graph.nodes[k].get("pr", 0))) for k in graph}
    #neg_log2_labels_pos = {k: (v[0] - .06, v[1]) for k, v in pos.items()}

    #bbox = dict(fc="black", ec="white", lw=2)

    #nx.draw_networkx_labels(
    #    graph, neg_log2_labels_pos, neg_log2_labels, bbox=bbox, font_color="white", font_size=11
    #)
    '''
    # Add labels to the graph for the Pr() values:
    precision_labels = {k: graph.nodes[k]["pr"] for k in graph}
    precision_labels_pos = {k: (v[0], v[1] + 0.04) for k, v in pos.items()}

    nx.draw_networkx_labels(
        graph, precision_labels_pos, precision_labels, font_color="green", font_size=10
    )
    # Add labels to the graph for the negative log2 values:
    neg_log2_labels = {k: -log2(graph.nodes[k]["pr"]) for k in graph}
    neg_log2_labels_pos = {k: (v[0], v[1] + 0.06) for k, v in pos.items()}
    nx.draw_networkx_labels(
        graph, neg_log2_labels_pos, neg_log2_labels, font_color="blue", font_size=11
    )
    '''
    # plt.show()
    return plt.figure(1)


def calculate_remaining_uncertainty(benchmark_graph: nx.DiGraph, prediction_graph: nx.DiGraph) -> float:
    """ Calculate the Remaining Uncertainty (RU) for a given prediction set.

    Parameters
    ----------
    benchmark_graph
        nx.DiGraph representing an annotation benchmark as binary data
        with proteins for indices and GO terms for columns
    prediction_dataframe
        nx.DiGraph representing annotation predictions as floating point
        numbers (probabilities) with proteins for indices and GO terms for columns

    Returns
    --------
    float

    Notes
    --------
    Remaining Uncertainty is simply the total information content of the nodes
    in the ontology that are contained in true annotation T, but not in the
    predicted annotation P.

    """
    benchmark_nodes = set(benchmark_graph.nodes)
    prediction_nodes = set(prediction_graph.nodes)
    ru_nodes = benchmark_nodes.difference(prediction_nodes)

    # print("RU NODES:")
    # for n in ru_nodes:
    #    print("\t", n, benchmark_graph.nodes[n].get("pr"))

    #info_content = sum(
    #    [-log2(benchmark_graph.nodes[node].get("pr")) for node in benchmark_graph if
    #     node in ru_nodes]
    #)
    info_content = sum(
        [benchmark_graph.nodes[node].get("pr") for node in benchmark_graph if
         node in ru_nodes]
    )
    return info_content


def calculate_misinformation(benchmark_graph: nx.DiGraph, prediction_graph: nx.DiGraph) -> float:
    """ Calculate the Misinformation (MI) for a given prediction set.

    Parameters
    ----------
    benchmark_graph
        networkx.DiGraph representing an annotation benchmark as binary data
        with proteins for indices and GO terms for columns
    prediction_graph
        networkx.DiGraph representing annotation predictions as floating point
        numbers (probabilities) with proteins for indices and GO terms for columns

    Returns
    --------
    float

    Notes
    ---------
    misinformation corresponds to the total information content of the nodes
    along incorrect paths in the prediction graph P.
    """

    prediction_nodes = set(prediction_graph.nodes)
    benchmark_nodes = set(benchmark_graph.nodes)
    mi_nodes = prediction_nodes.difference(benchmark_nodes)
    #print("MI NODES:")
    #for n in mi_nodes:
    #    print("\t", n, prediction_graph.nodes[n].get("pr"))

    #info_content = sum(
    #    [-log2(prediction_graph.nodes[node].get("pr")) for node in prediction_graph if
    #     node in mi_nodes and prediction_graph.nodes[node].get("pr") != 0]
    #)

    info_content = sum(
        [prediction_graph.nodes[node].get("pr") for node in prediction_graph if
         node in mi_nodes and prediction_graph.nodes[node].get("pr") != 0]
    )

    return info_content

