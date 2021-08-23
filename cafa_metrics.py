import pandas as pd

def get_remaining_uncertainity(
        predicted_terms: set, benchmark_terms: set, dag_ia_df: pd.DataFrame
) -> float:
    """ Calculates the RU (remaining uncertainity) for a given collection of predicted terms.
    RU is the sum of the information accretion values for all nodes in the benchmark NOT in the prediction
    the predicted_terms arg should be thresholded before calling this function.

    Per https://academic.oup.com/bioinformatics/article/29/13/i53/195366,
    we calculate RU(τ) (τ is a decision threshold) by:
    1. Sum the RU values for all proteins at the given threshold
    2. divide by the number of proteins in the benchmark
    """
    ru_terms = benchmark_terms - predicted_terms

    df = dag_ia_df.loc[ru_terms, 'ia']
    return df.sum()


    benchmark_diff_ia = []
    #test = [weighted_graph.nodes.get(term).get("ia", 0) for term in terms]

    benchmark_diff_ia = [
        weighted_graph.nodes[node].get("ia", 0)
        for node in weighted_graph.nodes
        if node in ru_terms
    ]

    return sum(benchmark_diff_ia)


def get_misinformation(
        predicted_terms: set, benchmark_terms: set, dag_ia_df: pd.DataFrame
) -> float:
    """ Calculates the MI (misinformation) for the given sets of predicted terms and benchmark terms.

    MI is the sum of the information accretion values for all nodes in the prediction NOT in the benchmark.

    The general use case is to pass the predicted and benchmark terms for an individual protein.
    """
    mi_terms = predicted_terms - benchmark_terms

    # mi_terms tends to be very large with many terms missing from the graph,
    # for efficiencies sake, get the subset that is actually in the graph:
    #all_graph_terms = set(weighted_graph.nodes.keys())

    #mi_ia_vals = [weighted_graph.nodes[term].get("ia") for term in (mi_terms & all_graph_terms)]
    #mi_ia_vals = [
    #    weighted_graph.nodes[node].get("ia", 0)
    #    for node in weighted_graph.nodes
    #    if node in mi_terms
    #]
    #dag_ia_df = pd.read_pickle("./data/weighted_dags/weighted_graph_CCO_dataframe.pkl")
    mi_terms_filtered = mi_terms & set(dag_ia_df.index)

    #print(dag_ia_df.loc[mi_terms_filtered, 'ia'].sum())

    df = dag_ia_df.loc[mi_terms_filtered, 'ia']
    return df.sum()
    #return 0 #sum(mi_ia_vals)


def get_rumi(
        predicted_terms: set, benchmark_terms: set, weighted_graph
) -> tuple:

    return (
        get_remaining_uncertainity(predicted_terms, benchmark_terms, weighted_graph),
        get_misinformation(predicted_terms, benchmark_terms, weighted_graph),
    )

