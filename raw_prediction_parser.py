from typing import Iterable, Optional
from datetime import datetime
import os
from pathlib import Path
import operator
import sys
import json
import numpy as np
import pandas as pd
from config import taxonomy_map, ontologies


def get_propagated_prediction_dataframe(prediction_df: pd.DataFrame, dag_df: pd.DataFrame):
    # Create a DataFrame joining the 'raw' prediction data with the DAG propagation DataFrame:
    term_ids = dag_df.columns.values
    prediction_df = pd.merge(prediction_df, dag_df, right_index=True, left_on='term', how='left')
    prediction_df.set_index('protein', inplace=True)

    # identify columns for propagation:
    propagation_mask = prediction_df[term_ids] == 1
    # Here we copy the threshold values to the propagated term columns and
    # simultaneously drop the term and threshold columns. There are still multiple
    # rows per protein:
    prediction_df = prediction_df[term_ids].mask(
        propagation_mask,
        prediction_df["threshold"],
        axis=0
    )

    # Aggregate the rows on protein, keeping the max threshold value for each protein/term pair:
    prediction_df = prediction_df.groupby("protein").aggregate("max")
    return prediction_df


def filter_dataframe(df: pd.DataFrame, filter_proteins:Optional[Iterable] = None, filter_terms: Optional[Iterable] = None):
    # filter the raw_prediction_df by the benchmark proteins and the benchmark ontology terms:

    if filter_proteins is not None:
        # TODO: Refactor the masking here to be less hokey:
        protein_mask = df['protein'].tolist()
        protein_mask = [v in filter_proteins for v in protein_mask]
        df = df[protein_mask]
    if filter_terms is not None:
        # TODO: This is fragile if the DataFrame columns change in anyway:
        terms_mask = [v in filter_terms for v in df.loc[:, 'term']]
        df = df[terms_mask]

    return df


def prediction_dataframe_to_dict(prediction_df: pd.DataFrame) -> dict:
    '''
    Converts a DataFrame of prediction data to a nested dictionary
    structure with this shape:
    {
       "T72270000115": {
          "GO:0000151": 0.01,
          "GO:0000228": 0.02,
          "GO:0000785": 0.02,
          "GO:0000790": 0.02,
          "GO:0005575": 0.3,
          ...
    '''

    prediction_dict = {}

    for protein, term_values in prediction_df.iterrows():
        mask = term_values > 0
        non_zero_annotations = term_values[mask]
        prediction_dict[protein] = non_zero_annotations.round(2).to_dict()

    return prediction_dict


if __name__ == "__main__":

    predictions_parent_directory = (
        "/media/scott/data/cafa3_submissions/ZhangFreddolinoLab/"
    )
    prediction_file_delimiter = "\t"
    benchmark_directory_path_str = "./data/benchmark/"
    benchmark_directory_path = Path(benchmark_directory_path_str)

    for ontology in ontologies:

        dag_df_filepath = f"../code/CLEAN/v6/data/propagation/propagation_map_df_{ontology.upper()}.pkl"
        dag_df = pd.read_pickle(dag_df_filepath)

        benchmark_files = list(benchmark_directory_path.glob(f"{ontology}_*.pkl"))

        print(f"PARSING {ontology}\n*********************\n")
        propagation_map_df_filepath = (
            f"../code/CLEAN/v6/data/propagation/propagation_map_df_{ontology}.pkl"
        )

        prediction_path = Path(predictions_parent_directory)
        prediction_files = prediction_path.glob("*.txt")

        for prediction_file in prediction_files:
            lab, model, taxon_id, *rest = prediction_file.stem.split("_")

            if int(taxon_id) not in taxonomy_map.keys():
                continue

            # find the relevant benchmark DataFrame:
            try:
                benchmark_file = list(benchmark_directory_path.glob(f"{ontology}_*_{taxon_id}_benchmark.pkl"))[0]
            except IndexError:
                # No benchmark found
                # TODO: Fix this exception, it's likely due to bad benchmark parsing
                continue

            taxon_ontology_benchmark_df = pd.read_pickle(benchmark_file)

            raw_prediction_df = pd.read_csv(
                prediction_file,
                engine="python",
                delimiter=prediction_file_delimiter,
                names=("protein", "term", "threshold"),
                skiprows=3,
                skipfooter=1,
                dtype={"protein": "string", "term": "string", "threshold": "float32"},
            )

            # TODO: This is fragile if the DataFrame columns change in anyway:
            benchmark_terms = taxon_ontology_benchmark_df.columns[2:]
            benchmark_proteins = taxon_ontology_benchmark_df.index.tolist()
            raw_prediction_df = filter_dataframe(
                raw_prediction_df,
                filter_proteins=benchmark_proteins,
                filter_terms=benchmark_terms
            )

            raw_prediction_df = get_propagated_prediction_dataframe(prediction_df=raw_prediction_df, dag_df=dag_df)

            '''
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
            '''
            # Convert the DataFrame to a more succint dictionary:
            prediction_dict = prediction_dataframe_to_dict(raw_prediction_df)

            output_filepath = f"./data/ZhangFreddolinoLab/{prediction_file.stem}_{ontology}.json"
            print(f"\tWRITING {output_filepath}")
            with open(output_filepath, 'w') as write_handle:
                json.dump(prediction_dict, write_handle)