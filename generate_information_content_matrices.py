""" Code for computing IA content values for graph nodes on a per-ontology basis.

This file generates pandas DataFrames with this general form:
+------------+----------+-------------+-------------+
|            |   weight |   precision |          ia |
+============+==========+=============+=============+
| GO:0005794 |       -5 |    0.365385 |  1.0068     |
+------------+----------+-------------+-------------+
| GO:0044431 |       -4 |    0.888889 |  0.117783   |
+------------+----------+-------------+-------------+
| GO:0043231 |       -4 |    1        | -0          |
+------------+----------+-------------+-------------+
| GO:0005622 |       -2 |    0.838164 |  0.176541   |
+------------+----------+-------------+-------------+
| GO:0044424 |       -2 |    0.991354 |  0.00868312 |
+------------+----------+-------------+-------------+
| GO:0031090 |       -2 |    0.44     |  0.820981   |
+------------+----------+-------------+-------------+
| GO:0044464 |       -1 |    1        | -0          |
+------------+----------+-------------+-------------+
| GO:0044446 |       -3 |    0.986755 |  0.0133335  |
+------------+----------+-------------+-------------+
| GO:0012505 |       -2 |    0.130435 |  2.03688    |
+------------+----------+-------------+-------------+
| GO:0016020 |       -1 |    0.29703  |  1.21392    |
+------------+----------+-------------+-------------+

"""
from pathlib import Path
from typing import Iterable
import pandas as pd
from goatools.obo_parser import GODag
import networkx as nx
from utils import parse_annotation_file, get_annotation_dataframe
from graph_utils import (
    get_graph,
    compute_nodes_information_content,
)


def get_ia_graph(
    obo_filepath: str, namespace_long: str, groundtruth_filepath: str, propagation_map_filepath: str
) -> nx.MultiDiGraph:
    """ Generates a DAG in the form of a networkx graph containing Information Content values for
    each graph node. """

    optional_attrs = ["relationship", "replaced_by", "consider"]
    namespaces = (namespace_long,)
    dag = GODag(
        obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None
    )

    # By using obj.item_id as the dictionary key instead of id,
    # we get canonical IDs instead of alt-IDs as keys and reduce
    # the size/redundancy of the dict:
    subdag = {
        obj.item_id: obj for id, obj in dag.items() if obj.namespace in namespaces
    }

    # sort the annotation pairs by protein ID only for the sake of consistency:
    benchmark_pairs = parse_annotation_file(groundtruth_filepath, sort=True)

    benchmark_df = get_annotation_dataframe(
        propagation_map_filepath=propagation_map_filepath,
        dag=subdag,
        annotation_pairs=benchmark_pairs,
        pare_df=True,
    )

    graph = get_graph(annotation_dataframe=benchmark_df, dag=subdag)
    compute_nodes_information_content(
        graph=graph, annotation_dataframe=benchmark_df, dag=subdag
    )

    return graph


def make_information_content_matrices(
    namespaces: Iterable,
    obo_filepath: str,
    output_path_str: str,
    verbose: str = True,
    output_filename_pattern=None,
) -> None:
    """Computes and writes pandas DataFrames representing the Information Content
    of a DAG on a per-ontology basis.

    The produced DataFrames have this general form:
    +------------+----------+-------------+-------------+
    |            |   weight |   precision |          ia |
    +============+==========+=============+=============+
    | GO:0005794 |       -5 |    0.365385 |  1.0068     |
    +------------+----------+-------------+-------------+
    | GO:0044431 |       -4 |    0.888889 |  0.117783   |
    +------------+----------+-------------+-------------+
    | ...        |      ... |         ... |       ...   |
    +------------+----------+-------------+-------------+
    """

    if output_filename_pattern is None:
        output_filename_pattern = "{namespace_short}_ia.pkl"

    output_path = Path(output_path_str)
    output_path.mkdir(parents=True, exist_ok=True)

    for namespace_dict in namespaces:
        (
            namespace_short,
            namespace_long,
            benchmark_filepath,
            propagation_map_filepath,
        ) = namespace_dict.values()

        if verbose:
            print(f"PROCESSING {namespace_short}/{namespace_long}")
            print(f"\tUSING benchmark {benchmark_filepath}")
            print(f"\tUSING {propagation_map_filepath}")

        weight_graph = get_ia_graph(
            obo_filepath, namespace_long, benchmark_filepath, propagation_map_filepath
        )
        graph_dict = {term: node for term, node in weight_graph.nodes.items()}
        dag_ia_df = pd.DataFrame.from_dict(data=graph_dict, orient="index")
        dag_id_df_filepath = output_path / output_filename_pattern.format(
            namespace_short=namespace_short
        )
        dag_ia_df.to_pickle(dag_id_df_filepath)

        if verbose:
            print(f"\tWRITING IA FOR {dag_ia_df.shape[0]} TERMS")
            print(f"\tWRITING TO {dag_id_df_filepath}")
            print("")


if __name__ == "__main__":
    import sys
    import yaml

    config_filepath = sys.argv[1]
    with open(config_filepath, "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        obo_filepath = config.get("obo_filepath")
        ontologies = config.get("ontologies")
        dag_directory = config.get("dag_directory")

        make_information_content_matrices(
            ontologies, obo_filepath, output_path_str=dag_directory
        )
