from typing import Iterable, Optional
from datetime import datetime
import os
from pathlib import Path
import operator
import sys
import json
import numpy as np
import pandas as pd
from config import taxonomy_map


def get_propagated_prediction_dataframe(
    prediction_df: pd.DataFrame, dag_df: pd.DataFrame
):
    # Create a DataFrame joining the 'raw' prediction data with the DAG propagation DataFrame:
    term_ids = dag_df.columns.values

    prediction_df = pd.merge(
        prediction_df, dag_df, right_index=True, left_on="term", how="left"
    )
    prediction_df.set_index("protein", inplace=True)

    # identify columns for propagation:
    propagation_mask = prediction_df[term_ids] == 1
    # Here we copy the threshold values to the propagated term columns and
    # simultaneously drop the term and threshold columns. There are still multiple
    # rows per protein:
    prediction_df = prediction_df[term_ids].mask(
        propagation_mask, prediction_df["threshold"], axis=0
    )

    # Aggregate the rows on protein, keeping the max threshold value for each protein/term pair:
    prediction_df = prediction_df.groupby("protein").aggregate("max")
    return prediction_df


def filter_dataframe(
    df: pd.DataFrame,
    filter_proteins: Optional[Iterable] = None,
    filter_terms: Optional[Iterable] = None,
):
    ''' FIlter the given DataFrame by the given proteins and the given ontology terms.
    It is assumed that the DataFrame has proteins has indices and GO terms as columns.

    '''

    if filter_proteins is not None:
        # TODO: Refactor the masking here to be less hokey:
        protein_mask = df["protein"].tolist()
        protein_mask = [v in filter_proteins for v in protein_mask]
        df = df[protein_mask]
    if filter_terms is not None:
        # TODO: This is fragile if the DataFrame columns change in anyway:
        terms_mask = [v in filter_terms for v in df.loc[:, "term"]]
        df = df[terms_mask]

    return df


def prediction_dataframe_to_dict(prediction_df: pd.DataFrame) -> dict:
    """
    Converts a DataFrame of prediction data (proteins as indices and GO terms as columns) to a nested dictionary
    structure with this general form:
    {
       "T72270000115": {
          "GO:0000151": 0.01,
          "GO:0000228": 0.02,
          "GO:0000785": 0.02,
          "GO:0000790": 0.02,
          "GO:0005575": 0.3,
          ...
    """

    prediction_dict = {}

    for protein, term_values in prediction_df.iterrows():
        mask = term_values > 0
        non_zero_annotations = term_values[mask]
        prediction_dict[protein] = non_zero_annotations.round(2).to_dict()

    return prediction_dict


def main(
    predictions_parent_directory: str,
    benchmark_parent_directory: str,
    dag_df_parent_directory: str,
    propagation_df_parent_directory: str,
    output_parent_directory: str,
    knowledge_type_id: int,
    ontologies: Iterable,
    prediction_file_delimiter: str,
    prediction_file_header_line_count: int = 3
):
    ''' Reads raw prediction data, filters the data based on benchmark data, propagates the annotations and writes
    the resulting filtered, propagated annotation predictions with probability thresholds to json-formatted files.
    '''


    predictions_parent_path = Path(predictions_parent_directory)
    benchmark_directory_path = Path(benchmark_parent_directory)
    dag_df_directory_path = Path(dag_df_parent_directory)
    output_directory_path = Path(output_parent_directory)
    propagation_directory_path = Path(propagation_df_parent_directory)
    written_files = []

    print(ontologies)

    for ontology in ontologies:

        dag_df_filepath = (
            propagation_directory_path / f"propagation_map_df_{ontology.upper()}.pkl"
        )
        dag_df = pd.read_pickle(dag_df_filepath)

        #benchmark_files = list(benchmark_directory_path.glob(f"{ontology}_*.pkl"))

        print(f"PARSING {ontology}\n*********************\n")

        prediction_files = list(predictions_parent_path.glob("*.txt"))
        prediction_files.sort(key=lambda f: f.stem)

        for prediction_file in prediction_files:
            print(f"\nPARSING {prediction_file}")
            lab, model, taxon_id, *rest = prediction_file.stem.split("_")

            # taxonomy_map is a dict mapping ncbi IDs (9606 for example) to uniprot shorthand (HUMAN for example).
            # benchmark files use the uniprot strings, predictions use the numerical IDs, so the mapping is necessary
            # for matching a prediction to a benchmark.
            if int(taxon_id) not in taxonomy_map.keys():
                print(f"\tSKIPPING. `{taxon_id}` NOT FOUND IN KNOWN TAXON IDS")
                continue

            # find the relevant benchmark DataFrame:
            try:
                benchmark_file = list(
                    benchmark_directory_path.glob(
                        f"{ontology}_*_{taxon_id}_type_{knowledge_type_id}_benchmark.json"
                    )
                )[0]
            except IndexError:
                # No benchmark found
                # TODO: Fix this exception, it's likely due to bad benchmark parsing
                print(f"\tNO BENCHMARK FOUND FOR {ontology} AND {taxon_id}")
                continue

            #taxon_ontology_benchmark_df = pd.read_pickle(benchmark_file)
            with open(benchmark_file, "r") as read_benchmark_handle:
                benchmark = json.load(read_benchmark_handle)

            # Not all CAFA 3 files seem to have the same header metadata, it's either 3 or 4 lines:
            print(f"\tSKIPPING {prediction_file_header_line_count} HEADER LINES")

            raw_prediction_df = pd.read_csv(
                prediction_file,
                engine="python",
                delimiter=prediction_file_delimiter,
                names=("protein", "term", "threshold"),
                skiprows=prediction_file_header_line_count,
                skipfooter=1,
                dtype={"protein": "string", "term": "string", "threshold": "float32"},
            )

            # TODO: This is fragile if the DataFrame columns change in anyway:
            #benchmark_terms = taxon_ontology_benchmark_df.columns[2:]
            #benchmark_proteins = taxon_ontology_benchmark_df.index.tolist()

            benchmark_terms = {v2 for v in benchmark.get("protein_annotations").values() for v2 in v}
            benchmark_proteins = benchmark.get("protein_annotations").keys()

            raw_prediction_df = filter_dataframe(
                raw_prediction_df,
                filter_proteins=benchmark_proteins,
                filter_terms=benchmark_terms,
            )

            # For empty predictions for the given benchmark, there is not reason to continue:
            if raw_prediction_df.shape[0] == 0:
                # This is an empty DataFrame, likely b/c the intersection of predicted proteins and benchmark proteins
                # is an empty set.
                print(f"\tSKIPPING EMPTY DATAFRAME FOR {ontology} AND {prediction_file}")
                continue

            raw_prediction_df = get_propagated_prediction_dataframe(
                prediction_df=raw_prediction_df, dag_df=dag_df
            )



            """
            At this point, we have a DataFrame with this form:
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            | protein       |   GO:0000015 |   GO:0000109 |   GO:0000110 |   GO:0000111 |   GO:0000112 |   GO:0000113 |
            +===============+==============+==============+==============+==============+==============+==============+
            | T100900000026 |            0 |            0 |            0 |            0 |            0 |            0 |
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            | T100900000115 |            0 |            0 |            0 |            0 |            0 |            0 |
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            | T100900000116 |            0 |            0 |            0 |            0 |            0 |            0 |
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            | T100900000167 |            0 |            0 |            0 |            0 |            0 |            0 |
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            | T100900000453 |            0 |            0 |            0 |            0 |            0 |            0 |
            +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
            """
            # Convert the DataFrame to a more succinct dictionary:
            prediction_dict = prediction_dataframe_to_dict(raw_prediction_df)

            # output_filepath = f"./data/ZhangFreddolinoLab/{prediction_file.stem}_{ontology}.json"
            #output_directory_path = Path(output_directory)
            output_directory_path.mkdir(parents=True, exist_ok=True)
            #output_filepath = output_directory_path / f"{prediction_file.stem}_{ontology}.json"
            output_filepath = (
                output_directory_path / f"{prediction_file.stem}_{ontology}.json"
            )
            print(f"\tWRITING {output_filepath}")
            written_files.append(output_filepath)
            with open(output_filepath, "w") as write_handle:
                json.dump(prediction_dict, write_handle)

        print(f"FINISHED PARSING {ontology}\n*********************\n")
    return written_files

if __name__ == "__main__":
    '''
    predictions_parent_directory = "/media/airport/cafa3_submissions/ZhangFreddolinoLab"
    prediction_file_delimiter = "\t"
    benchmark_directory_path_str = "./data/parsed_benchmark"
    # benchmark_directory_path = Path(benchmark_directory_path_str)
    dag_directory = "./data/dag_ia"
    # main(predictions_parent_directory, benchmark_directory_path, )
    knowledge_type = 1
    '''

    import yaml
    import sys

    config_filepath = sys.argv[1]
    print(f"USING CONFIG {config_filepath}")

    with open(config_filepath, "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        ontologies = config.get("ontologies")
        benchmark_json_directory = config.get("benchmark_json_directory")
        dag_ic_directory = config.get("dag_ic_directory")
        prediction_file_delimiter = config.get("prediction_file_delimiter", " ")
        knowledge_type = config.get("knowledge_type")
        propagation_df_directory = config.get("propagation_df_directory")
        skip_lines = int(config.get("prediction_file_header_line_count", 3))

        raw_predictions_directory = config.get("raw_predictions_directory")
        output_directory = config.get("predictions_json_directory")

        parser_output = main(
            raw_predictions_directory,
            benchmark_json_directory,
            dag_ic_directory,
            propagation_df_directory,
            output_directory,
            knowledge_type,
            ontologies,
            prediction_file_delimiter,
            skip_lines
        )

        print(f"\nCREATED FILES:")
        for f in parser_output:
            print(f"\t{f}")
