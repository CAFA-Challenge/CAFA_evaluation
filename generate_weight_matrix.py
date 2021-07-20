""" This file is a test/proof-of-concept for
generating graph nodes information content
weights.
"""
import sys
from math import log, isclose
from typing import Iterable
import pandas as pd
from goatools.obo_parser import GODag
import networkx as nx
from utils import parse_annotation_file, get_annotation_dataframe
from graph_utils import (
    get_graph,
    compute_nodes_information_content,
)

'''
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


def calculate_node_precision(
    term_id: str, dag: GODag, annotation_dataframe: pd.DataFrame
) -> float:
    """ Calculates the 'precision' value for a graph node (annotation) based on the frequency
    of the annotation in a dataset compared to the frequency of annotation for the node's
    direct parent nodes.

    """
    parent_terms = get_parent_terms(
        term_id=term_id, dag=dag, relationships=("part_of",)
    )
    # this is a boolean matrix of all proteins x the GO terms of interest
    # (terms that are parents of the term being evaluated:
    parent_mask = annotation_dataframe.loc[:, parent_terms].all(axis=1)
    parent_count = len(annotation_dataframe.loc[parent_mask]) + 1
    term_count = sum(annotation_dataframe.loc[:, term_id]) + 1
    return round(-log(term_count / parent_count), 2)
'''

def main(namespace_short, namespace_long):
    # matlab_df = pd.read_pickle("matlab_cco_node_weight_matrix.pkl")
    #matlab_df = pd.read_csv("../v3/data/MATLAB/matlab_working_weights_cco.txt")
    #matlab_df.set_index("go_term", inplace=True, drop=True)

    optional_attrs = ["relationship", "replaced_by", "consider"]
    obo_filepath = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    namespaces = (namespace_long,)
    # predictions_df_path = "./kiharalab1_1_9606_working.pkl"
    #predictions_df_path = "./data/predictions/ZZZ_1_9606_CCO.pkl"
    #protein_of_interest = "T96060014271"

    dag = GODag(
        obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None
    )

    ontology = namespace_short #"MFO"
    groundtruth_filepath = f"/home/scott/virtualenvs/cafa_scott/code/CAFA_assessment_tool/data/benchmark/CAFA3_benchmarks/groundtruth/leafonly_{ontology}.txt"

    propagation_map_filepath = (
        f"../v3/data/propagation/propagation_matrix_dataframe_{ontology}_FULL.pkl"
    )

    #propagation_map_filepath = (
    #    f"./data/propagation/propagation_matrix_dataframe_{ontology}_is_a_part_of.pkl"
    #)
    prop_map = pd.read_pickle(propagation_map_filepath)

    # propagation_map_filepath = "./data/cco_propagation_map_is_a_only.pkl"

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
    compute_nodes_information_content(graph=graph, annotation_dataframe=benchmark_df, dag=subdag)

    # pare the benchmark down to the relevant GO terms (columns):
    terms_of_interest = list(graph.nodes)
    #benchmark_df = benchmark_df.loc[:, terms_of_interest]

    matlab_data = []
    # read in test data from matlab:
    matlab_filepath = f"/home/scott/Documents/MATLAB/CAFA2/matlab_support_sums_{ontology.lower()}.txt"
    with open(matlab_filepath, "r") as matlab_read_handle:
        # ignore the header line:
        matlab_read_handle.readline()
        for line in matlab_read_handle:
            line_split = line.rstrip().split(",")
            matlab_data.append(
                (
                    line_split[1],
                    int(line_split[2]),
                    int(line_split[-4]),
                    float(line_split[-2]),
                    float(line_split[-1]),
                )
            )

        matlab_read_handle.seek(0)
        matlab_read_handle.seek(0)

    passed_tests = 0
    failed_tests = 0
    missing_go_terms = 0


    # Iterate over the Matlab GO Terms testing if each term in the Matlab data exists in the networkx graph:
    for line in matlab_data:
        ml_term_id, ml_parent_count, ml_term_count, ml_precision, ml_nlp = line
        assert ml_term_id in graph.nodes
        continue

        info_content = float(graph.nodes[ml_term_id].get("ia"))
        print(ml_term_id, "ML:", ml_nlp, "PY:", info_content, ml_nlp == info_content)
        print([graph.nodes[node] for node in graph.nodes if node == ml_term_id])
        parents = [list(graph.predecessors(node)) for node in graph.nodes if node == ml_term_id]
        # flatten the nested list(s):
        parents = [item for sublist in parents for item in sublist]
        print(parents)
        print([graph.nodes[parent] for parent in parents])
        #assert ml_nlp == info_content
        print("===================")

    matlab_terms = set([row[0] for row in matlab_data])
    print(matlab_terms)

    # test the count of proteins annotated with each node (term) in the graph
    for term in graph.nodes:
        try:
            ml_term_id, ml_parent_count, ml_term_count, ml_precision, ml_nlp = [
                row for row in matlab_data if row[0] == term
            ][0]



        except IndexError:
            # go term not found in Matlab data:
            print(f"{term} NOT found in Matlab data")
            print(f"\t{graph.nodes[term]}")

            prop_row = prop_map.loc[term,:]
            mask = prop_row == 1
            prop_row = prop_row.loc[mask]
            for k in prop_row.keys():
                print("\t", k, k in matlab_terms)

            missing_go_terms += 1
            continue

        info_content = float(graph.nodes[term].get("ia"))
        #info_content = calculate_node_precision(
        #    term_id=term, dag=dag, annotation_dataframe=benchmark_df
        #)

        try:
            assert isclose(info_content, ml_nlp, abs_tol=0.001)
            passed_tests += 1
        except AssertionError:
            print(term, info_content, ml_nlp)
            failed_tests += 1

    print(f"{missing_go_terms:,} TESTS SKIPPED, {passed_tests:,} TESTS PASSED, {failed_tests:,} TESTS FAILED")

    return graph


if __name__ == "__main__":
    namespace_short = "CCO"
    namespace_long = "cellular_component"
    weight_graph = main(namespace_short, namespace_long)
    nx.write_gpickle(weight_graph, f"./data/weighted_graph_{namespace_short}.pkl")

    for term, node in weight_graph.nodes.items():
        print(term, node)
