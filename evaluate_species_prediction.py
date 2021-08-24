from pathlib import Path
from typing import Iterable
import json
import numpy as np
import pandas as pd
import yaml
from cafa_metrics import get_rumi


def initialize_proteins_and_thresholds_dataframe(
    proteins: Iterable, thresholds: Iterable
) -> pd.DataFrame:
    """Creates a pandas.DataFrame that will have one row per protein/threshold pair.
    # It will only include the thresholds that actually exist in the prediction data and not a
    # complete range of possible thresholds:
    :param proteins:
    :param thresholds:
    :return:
    """
    columns = (
        "protein",
        "threshold",
        "tp_ia",
        "fp_ia",
        "fn_ia",
        "benchmark_ia",
        "weighted_precision",
        "weighted_recall",
        "tp",
        "fp",
        "fn",
        "tn",
        "precision",
        "recall",
    )

    matrix = np.zeros(((len(proteins) * len(thresholds)), len(columns)))
    protein_and_threshold_df = pd.DataFrame(
        data=matrix,
        columns=columns
    )
    df_datatypes = {
        "protein": "str",
        "threshold": "float",
        "tp_ia": "float",
        "fp_ia": "float",
        "fn_ia": "float",
        "benchmark_ia": "float",
        "weighted_precision": "float",
        "weighted_recall": "float",
        "tp": "float",
        "fp": "float",
        "fn": "float",
        "tn": "float",
        "precision": "float",
        "recall": "float",
    }
    protein_and_threshold_df = protein_and_threshold_df.astype(df_datatypes)

    # The DataFrame is empty. Next, populate the protein and threshold columns:
    benchmark_protein_index = -1

    for i in range(0, protein_and_threshold_df.shape[0]):
        threshold_index = i % len(thresholds)

        if threshold_index == 0:
            # for every complete loop of the threshold range,
            # increment the protein index:
            benchmark_protein_index += 1

        protein_and_threshold_df.iloc[i, 0] = proteins[benchmark_protein_index]
        protein_and_threshold_df.iloc[i, 1] = thresholds[threshold_index]

    protein_and_threshold_df.set_index(
        ["protein", "threshold"], drop=True, inplace=True
    )

    """ At this point, protein_and_threshold_df has this form:
    +------------------------+------+------+------+------+
    |                        |   tp |   fp |   fn |   tn |
    +========================+======+======+======+======+
    | ('T72270000115', 0.01) |    0 |    0 |    0 |    0 |
    +------------------------+------+------+------+------+
    | ('T72270000115', 0.02) |    0 |    0 |    0 |    0 |
    +------------------------+------+------+------+------+
    | ('T72270000115', 0.03) |    0 |    0 |    0 |    0 |
    +------------------------+------+------+------+------+
    | ('T72270000115', 0.05) |    0 |    0 |    0 |    0 |
    +------------------------+------+------+------+------+
    """
    return protein_and_threshold_df


def get_confusion_matrix_terms(predicted_terms: set, benchmark_terms: set) -> dict:
    """For two given sets of terms, compute and return:
    * terms for the true positive set
    * terms for the false positive set
    * terms for the false_negative set
    """
    true_positive_terms = predicted_terms & benchmark_terms
    false_positive_terms = predicted_terms - benchmark_terms
    false_negative_terms = benchmark_terms - predicted_terms

    return {
        "TP": true_positive_terms,
        "FP": false_positive_terms,
        "FN": false_negative_terms,
    }


def calculate_weighted_confusion_matrix(
    predicted_terms: set, benchmark_terms: set, ia_df: pd.DataFrame
) -> dict:
    """Calculates the weighted precision and recall for two sets of terms
    Weighted precision and recall rely on the information content (IC) of relevant terms (nodes).
    Here we retrieve the IC for the relevant nodes from the node_weights_df.

    In practice, this will be called for the predicted_terms @ a specific threshold.
    """
    cm_terms = get_confusion_matrix_terms(predicted_terms, benchmark_terms)
    tp_info_accretion = sum([ia_df.loc[term, "ia"] for term in cm_terms.get("TP", [])])
    fp_info_accretion = sum([ia_df.loc[term, "ia"] for term in cm_terms.get("FP", [])])
    fn_info_accretion = sum([ia_df.loc[term, "ia"] for term in cm_terms.get("FN", [])])
    benchmark_info_accretion = ia_df.loc[benchmark_terms, "ia"].sum()

    try:
        weighted_precision = tp_info_accretion / (tp_info_accretion + fp_info_accretion)
    except ZeroDivisionError:
        weighted_precision = 0

    weighted_recall = tp_info_accretion / benchmark_info_accretion

    #print(get_rumi(predicted_terms, benchmark_terms, ia_df))
    #predicted_terms: set, benchmark_terms: set, weighted_graph

    return {
        "TP": tp_info_accretion,
        "FP": fp_info_accretion,
        "FN": fn_info_accretion,
        "benchmark_ia": benchmark_info_accretion,
        "weighted_precision": weighted_precision,
        "weighted_recall": weighted_recall
    }


def calculate_confusion_matrix(predicted_terms: set, benchmark_terms: set) -> dict:
    """Calculates true positive, false positive and false negative for two sets of terms.
    Does not calculate true negative.

    In practice, this will be called for the predicted_terms @ a specific threshold.
    """
    cm_terms = get_confusion_matrix_terms(predicted_terms, benchmark_terms)
    true_positive = len(cm_terms["TP"])
    false_positive = len(cm_terms["FP"])
    false_negative = len(cm_terms["FN"])
    try:
        precision = true_positive / (true_positive + false_positive)
    except ZeroDivisionError:
        precision = 0

    return {
        "TP": true_positive,
        "FP": false_positive,
        "FN": false_negative,
        "precision": precision,
        "recall": true_positive / len(benchmark_terms),
    }


def get_confusion_matrix_dataframe(
    prediction_dict: dict, benchmark_dict: dict, ia_df: pd.DataFrame
) -> pd.DataFrame:
    """Constructs a pandas.DataFrame with a row for each protein/threshold pair.
    The proteins are sourced from the benchmark_dict and the thresholds are sourced
    from the prediction dict.

    The prediction_dict maps proteins to terms and threshold values and should have this form:
    {
        "T72270000115": {
            "GO:0000151": 0.01, "GO:0000228": 0.02, "GO:0000785": 0.02, "GO:0000790": 0.02, "GO:0005575": 0.3, ...
        },
        "T72270000700": {
            "GO:0000151": 0.01, "GO:0000307": 0.02, "GO:0000428": 0.02, "GO:0005575": 0.07, "GO:0005576": 0.01, ...
        },

    And the benchmark_dict should have this form:
    {
        "benchmark_taxon": "DROME",
        "benchmark_taxon_id": null,
        "benchmark_ontology": "bpo",
        "benchmark_ontology_term_count": 28678,
        "protein_annotations": {
            "T72270000015": ["GO:0007154", "GO:0007165", "GO:0007186", "GO:0007602", "GO:0007603", ...],
             "T72270000115": ["GO:0000003", "GO:0007276", "GO:0007283", "GO:0008150", "GO:0010468", ...],
        }
    }

    """

    # Get all threshold values from the nested dictionaries in the prediction data:
    distinct_prediction_thresholds = sorted(
        list(
            {
                threshold
                for protein in prediction_dict.values()
                for threshold in protein.values()
            }
        )
    )

    # the benchmark json data should include this metadata:
    benchmark_ontology = benchmark_dict.get("benchmark_ontology")
    benchmark_ontology_term_count = benchmark_dict.get("benchmark_ontology_term_count")
    benchmark_taxon = benchmark_dict.get("benchmark_taxon")
    benchmark_taxon_id = benchmark_dict.get("benchmark_taxon_id")
    benchmark_annotations = benchmark_dict.get("protein_annotations")
    benchmark_proteins = list(benchmark_annotations.keys())

    # Next, create a pandas.DataFrame that will have one row per protein/threshold pair.
    # It will only include the thresholds that actually exist in the prediction data and not a
    # complete range of possible thresholds:

    protein_and_threshold_df = initialize_proteins_and_thresholds_dataframe(
        proteins=benchmark_proteins, thresholds=distinct_prediction_thresholds
    )

    # protein_and_threshold_df has keys (proteins and thresholds), but no values.
    # Next, populate the DataFrame with the confusion matrix values @ each
    # decision threshold:
    for threshold in distinct_prediction_thresholds:
        for protein in benchmark_proteins:

            predicted_terms = prediction_dict.get(protein, {})
            # Limit the predictions by the threshold at hand:
            predicted_annotations = {
                k for k, v in predicted_terms.items() if v >= threshold
            }
            benchmark_protein_annotation = set(benchmark_annotations.get(protein))

            confusion_matrix = calculate_confusion_matrix(
                predicted_terms=predicted_annotations,
                benchmark_terms=benchmark_protein_annotation,
            )
            # print(confusion_matrix)
            protein_and_threshold_df.loc[protein, threshold].tp = confusion_matrix["TP"]
            protein_and_threshold_df.loc[protein, threshold].fp = confusion_matrix["FP"]
            protein_and_threshold_df.loc[protein, threshold].fn = confusion_matrix["FN"]
            true_negative = benchmark_ontology_term_count - sum(
                confusion_matrix.values()
            )
            protein_and_threshold_df.loc[protein, threshold].tn = true_negative

            protein_and_threshold_df.loc[
                protein, threshold
            ].precision = confusion_matrix["precision"]
            protein_and_threshold_df.loc[protein, threshold].recall = confusion_matrix[
                "recall"
            ]

            #print(threshold, protein)
            ia_sums = calculate_weighted_confusion_matrix(
                predicted_terms=predicted_annotations,
                benchmark_terms=benchmark_protein_annotation,
                ia_df=ia_df,
            )
            protein_and_threshold_df.loc[protein, threshold].tp_ia = ia_sums["TP"]
            protein_and_threshold_df.loc[protein, threshold].fp_ia = ia_sums["FP"]
            protein_and_threshold_df.loc[protein, threshold].fn_ia = ia_sums["FN"]
            protein_and_threshold_df.loc[protein, threshold].benchmark_ia = ia_sums[
                "benchmark_ia"
            ]
            protein_and_threshold_df.loc[
                protein, threshold
            ].weighted_precision = ia_sums["weighted_precision"]
            protein_and_threshold_df.loc[
                protein, threshold
            ].weighted_recall = ia_sums["weighted_recall"]

    # Lastly, add some metadata to each row:
    protein_and_threshold_df.insert(0, "taxon", benchmark_taxon)
    protein_and_threshold_df.insert(0, "taxon_id", benchmark_taxon_id)
    protein_and_threshold_df.insert(0, "ontology", benchmark_ontology)

    """ The final DataFrame has this form:
    +------------------------+------------+------------+---------+------+------+------+------+
    |                        | ontology   |   taxon_id | taxon   |   tp |   fp |   fn |   tn |
    +========================+============+============+=========+======+======+======+======+
    | ('T72270000115', 0.01) | CCO        |       7227 | DROME   |   10 |  102 |    0 | 3793 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.02) | CCO        |       7227 | DROME   |   10 |   44 |    0 | 3851 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.03) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.05) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.06) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+
    """
    return protein_and_threshold_df


def evaluate_species(
    prediction_filepath_str: str, benchmark_filepath_str: str, ia_df: pd.DataFrame
) -> pd.DataFrame:
    with open(prediction_filepath_str, "r") as prediction_handle:
        predictions = json.load(prediction_handle)

    with open(benchmark_filepath_str, "r") as benchmark_handle:
        benchmark = json.load(benchmark_handle)

    metrics_df = get_confusion_matrix_dataframe(
        prediction_dict=predictions, benchmark_dict=benchmark, ia_df=ia_df
    )
    return metrics_df


def main(
    predictions_parent_directory: str,
    benchmark_parent_directory: str,
    ia_parent_directory: str,
    model_id: int,
    ontologies: Iterable = ("CCO", "BPO"),
):
    """A generator function which yields pandas DataFrames. Each DataFrame contains evaluation metrics
    for a specific species + ontology pairing.

    predictions_parent_directory is a str representation of the directory which contains the 
    json-representation of the prediction data. These json files should have names derived from the raw prediction files. 
    Here's an example filename: DoeLab_1_9606_go_CCO.json
    That filename contains the original lab name (DoeLab), model number (1), a species ID (9606) and ontology shorthand (CCO). 

    For the prediction files, the file contents are nested dictionaries mapping proteins to terms with threshold values:

    {
        "T96060000375": {
            "GO:0000139": 0.02, 
            "GO:0000151": 0.01, 
            "GO:0005575": 0.06, 
            "GO:0005576": 0.02, 
            ...
        },
        "T96060000845": {
            "GO:0000151": .02,
            ...
        }
    }

    The yielded DataFrames have the following form with one row per protein + threshold pairing:
    +------------------------+------------+------------+---------+------+------+------+------+
    |                        | ontology   |   taxon_id | taxon   |   tp |   fp |   fn |   tn |
    +========================+============+============+=========+======+======+======+======+
    | ('T72270000115', 0.01) | CCO        |       7227 | DROME   |   10 |  102 |    0 | 3793 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.02) | CCO        |       7227 | DROME   |   10 |   44 |    0 | 3851 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.03) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.05) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+
    | ('T72270000115', 0.06) | CCO        |       7227 | DROME   |   10 |    3 |    0 | 3892 |
    +------------------------+------------+------------+---------+------+------+------+------+

    """

    predictions_path = Path(predictions_parent_directory)
    benchmark_path = Path(benchmark_parent_directory)
    ia_directory_path = Path(ia_parent_directory)

    for ontology in ontologies:
        ia_df = pd.read_pickle(ia_directory_path / f"{ontology}_ia.pkl")
        prediction_files = list(predictions_path.glob(f"*{model_id}_*{ontology}*json"))
        benchmark_files = list(benchmark_path.glob(f"*{ontology}*json"))

        for benchmark_file in benchmark_files:
            _, taxon_str, taxon_id, *rest = benchmark_file.stem.split("_")
            try:
                taxon_id = int(taxon_id)
            except ValueError:
                # No taxon_id found in file name
                continue

            prediction_file = [f for f in prediction_files if str(taxon_id) in str(f)]
            # should only find one file:
            if len(prediction_file) == 1:
                prediction_file = prediction_file[0]
            else:
                continue

            yield evaluate_species(prediction_file, benchmark_file, ia_df)


if __name__ == "__main__":

    with open("./parser_config.yml", "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        #prediction_filepath_str = config.get("prediction_filepath")
        prediction_filepath_str = config.get("predictions_json_directory")

        benchmark_filepath_str = config.get("benchmark_filepath")
        dag_directory_filepath_str = config.get("dag_directory")
        model_id = config.get("model_id")
        ontologies = config.get("ontologies")
        predictor_group_name = config.get("predictor_group_name")

        for taxon_result_df in main(
            prediction_filepath_str,
            benchmark_filepath_str,
            dag_directory_filepath_str,
            model_id,
            ontologies,
        ):
            print(taxon_result_df.iloc[:10, :].to_markdown(tablefmt="grid"))
            taxon = taxon_result_df.iloc[0, :].taxon
            ontology = taxon_result_df.iloc[0, :].ontology

            output_directory = Path(f"./data2/working/{predictor_group_name}")
            output_directory.mkdir(parents=True, exist_ok=True)

            taxon_result_df.to_pickle(
                output_directory / f"{taxon}_{ontology}_{model_id}.pkl"
            )

            print("\n\n")
