from typing import Iterable
from pathlib import Path
import numpy as np
import pandas as pd


def main(input_files: Iterable) -> pd.DataFrame:
    """Computes evaluation average metrics across species

    Generates a pandas DataFrame with this form:
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |   threshold | ontology   | taxons            |   average_precision |   average_recall |   average_weighted_precision |   average_weighted_recall |
    +=============+============+===================+=====================+==================+==============================+===========================+
    |        0.01 | cco        | HUMAN,MOUSE,ARATH |           0.270375  |        0.207061  |                    0.165135  |                0.116404   |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        0.02 | cco        | HUMAN,MOUSE,ARATH |           0.271606  |        0.1975    |                    0.166157  |                0.105961   |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        ...  | ...        | ...               |           ...       |        ...       |                    ...       |                ...        |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        0.97 | cco        | HUMAN,MOUSE,ARATH |           0.025446  |        0.0193828 |                    0.0161426 |                0.0100999  |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        0.98 | cco        | MOUSE,ARATH       |           0.0220311 |        0.0178682 |                    0.0141127 |                0.00957173 |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        0.99 | cco        | MOUSE,ARATH       |           0.0204188 |        0.0170323 |                    0.0136922 |                0.00928143 |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+
    |        1    | cco        | ARATH             |           0.0185571 |        0.0164021 |                    0.0145303 |                0.00983286 |
    +-------------+------------+-------------------+---------------------+------------------+------------------------------+---------------------------+


    """
    # Collect all thresholds across all of the species predictions:
    all_thresholds = set()

    for file in input_files:
        df = pd.read_pickle(file)
        thresholds = df.index.get_level_values(1)
        all_thresholds.update(thresholds)

    all_thresholds = sorted(list(all_thresholds))
    # At this point all_thresholds is a sorted list of unique prediction thresholds

    # Construct a dataframe to hold the computed results:
    columns = (
        "threshold",
        "ontology",
        #"taxons",
        "average_precision",
        "average_recall",
        "average_weighted_precision",
        "average_weighted_recall",
    )
    matrix = np.zeros((len(all_thresholds), len(columns)))
    results_df = pd.DataFrame(data=matrix, columns=columns)
    results_df = results_df.astype({"ontology": "str"})
    results_df["threshold"] = all_thresholds
    results_df.set_index("threshold", drop=True, inplace=True)

    for threshold in all_thresholds:
        all_species_at_threshold_df = None

        """ Read the species-specific DataFrames generated by evaluate_species_prediction.py and
        concatenate the data from each species on a per-threshold basis.
        """
        for file in input_files:

            df = pd.read_pickle(file)

            threshold_mask = df.index.get_level_values(1) >= threshold
            df = df[threshold_mask]
            if all_species_at_threshold_df is None:
                all_species_at_threshold_df = df
            else:
                all_species_at_threshold_df = all_species_at_threshold_df.append(df)


        protein_count = all_species_at_threshold_df.shape[0]

        average_precision = (
            all_species_at_threshold_df.loc[:, "precision"].sum() / protein_count
        )
        average_recall = (
            all_species_at_threshold_df.loc[:, "recall"].sum() / protein_count
        )
        average_weighted_precision = (
            all_species_at_threshold_df.loc[:, "weighted_precision"].sum()
            / protein_count
        )
        average_weighted_recall = (
            all_species_at_threshold_df.loc[:, "weighted_recall"].sum() / protein_count
        )

        average_ru = (
                all_species_at_threshold_df.loc[:, "ru"].sum() / protein_count
        )
        average_mi = (
                all_species_at_threshold_df.loc[:, "mi"].sum() / protein_count
        )

        results_df.loc[threshold, "ontology"] = all_species_at_threshold_df.iloc[
            0
        ].ontology
        results_df.loc[threshold, "average_precision"] = average_precision
        results_df.loc[threshold, "average_recall"] = average_recall
        results_df.loc[
            threshold, "average_weighted_precision"
        ] = average_weighted_precision
        results_df.loc[threshold, "average_weighted_recall"] = average_weighted_recall
        results_df.loc[threshold, "average_ru"] = average_ru
        results_df.loc[threshold, "average_mi"] = average_mi

        #taxons = ",".join(set(all_species_at_threshold_df.loc[:, "taxon"]))
        #results_df.loc[threshold, "taxons"] = taxons

    return results_df


if __name__ == "__main__":

    ''' This is dependent on the species/ontology specific DataFrames created by evaluate_species_prediction.py '''

    import sys
    import yaml
    config_filepath = sys.argv[1]

    with open(config_filepath, "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        dataframe_read_directory = Path(config.get("predictions_dataframes_directory"))
        ontologies = config.get("ontologies")

    for ontology in ontologies:
        input_filename_pattern = f"*_{ontology.lower()}_*.pkl"
        species_specific_files = list(dataframe_read_directory.glob(input_filename_pattern))
        print(f"\n{ontology}")
        print("READING")
        for ssf in species_specific_files:
            print(f"\t{ssf}")

        metrics_df = main(input_files=species_specific_files)
        print(metrics_df.to_markdown(tablefmt='grid'))
        print("\n")

