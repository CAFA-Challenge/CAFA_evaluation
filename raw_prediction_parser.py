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


def filter_prediction_dataframe(
    prediction_df: pd.DataFrame,
    proteins: Optional[Iterable] = None,
    terms: Optional[Iterable] = None,
):
    """ Filters a pandas.DataFrame with the assumed shape:
     0   protein    745 non-null    string
     1   term       745 non-null    string
     2   threshold  745 non-null    float32
    """
    if terms is not None:
        # limit the GO terms by ontology:
        prediction_df = prediction_df.loc[prediction_df["term"].isin(terms)]

    if proteins is not None:
        # limit the proteins by the benchmark:
        prediction_df = prediction_df.loc[prediction_df["protein"].isin(proteins)]

    return prediction_df


def get_propagated_dataframe(prediction_df, benchmark_proteins, dag_df):
    propagated = np.zeros(
        shape=(len(benchmark_proteins), len(dag_df.index)), dtype="ubyte"
    )

    # sort by threshold value so high threshold predictions are parsed last
    # in order to avoid overwriting high-probability predictions with low-probability
    # predictions
    prediction_df.sort_values(by=["threshold"], inplace=True)

    for _, protein, term, threshold in prediction_df.itertuples():
        # get the propagation terms from the dag:
        prop_mask = dag_df.loc[term, :] == 1
        # and the zero-based indices of the terms from the mask:
        propagation_indices = dag_df.columns.get_indexer(dag_df.columns[prop_mask])
        # get the zero-based index from the benchmark DataFrame of the protein:
        protein_index = benchmark_proteins.index(protein)
        # update the numpy array with the threshold data, converting the threshold values
        # to ints from zero to 1 in order to use a more efficient datatype.
        # because the incoming prediction data is sorted by threshold in ascending order,
        # we do not take care to overwrite higher values:
        propagated[protein_index, propagation_indices] = round(threshold, 2) * 100

    propagated = pd.DataFrame(
        data=propagated, index=benchmark_proteins, columns=dag_df.columns
    )
    propagated = propagated.astype(pd.SparseDtype("float", 0))
    propagated = propagated / 100
    return propagated


def parse_raw_predictions(
    predictions_parent_directory: str,
    benchmark_df_filepath: str,
    propagation_df_filepath: str,
    ontology_name: str,
    output_directory: str = "./",
    prediction_file_delimiter: str = " ",
    species_id_mapping: Optional[dict] = None
):
    ''' Parses the raw species-specific csv files in the given predictions_parent_directory into pickled
    pandas DataFrames. There will be one output pickle file per raw input file.

    :param predictions_parent_directory:
    :param benchmark_df_filepath:
    :param propagation_df_filepath:
    :param ontology_name:
    :param output_directory:
    :param prediction_file_delimiter:
    :param species_id_mapping:
    :return:
    '''

    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    benchmark_df = pd.read_pickle(benchmark_df_filepath)
    benchmark_proteins = benchmark_df.index
    del benchmark_df

    ontology_dag_df = pd.read_pickle(propagation_df_filepath)
    ontology_terms = ontology_dag_df.index
    # del ontology_dag_df

    # Collect the file object references:
    predictions_root_directory = Path(predictions_parent_directory)
    model_id = 1
    prediction_files = predictions_root_directory.glob(f"*_{model_id}_*.txt")

    # make a generator for tuples of file path and size: ('/Path/to/the.file', 1024)
    files_and_sizes = list(((path, os.path.getsize(path)) for path in prediction_files))

    sorted_files_with_size = sorted(
        files_and_sizes, key=operator.itemgetter(1), reverse=True
    )

    for f, s in sorted_files_with_size:
        start_t = datetime.now()
        print(f"\nPROCESSING {f}")
        lab, model_id, species_id, *rest = f.stem.split("_")
        species_id = int(species_id)

        #species_name = species_id_mapping.get(int(species))

        #print(species, species_name)


        knowledge_type_id = 1
        ontology = 'CCO'
        # Get the relevant proteins for the ontology and evaluation type:
        foobar_filepath = f"./data/{ontology}_proteins_by_species_and_knowledge_type.pkl"
        relevant_proteins_df = pd.read_pickle(foobar_filepath)

        knowledge_type_mask = relevant_proteins_df.loc[:, 'knowledge_type'] == knowledge_type_id
        relevant_proteins_df[knowledge_type_mask]
        species_id_mask = relevant_proteins_df.loc[:, 'species_id'] == species_id
        relevant_proteins_df = relevant_proteins_df[species_id_mask]
        if relevant_proteins_df.empty:
            print("\tSKIPPING PROCESS FOR FILE. NO SPECIES MATCH")
            print(f"\t{lab} {model_id}, {species_id}\n")
            continue


        relevant_benchmark_proteins = list(relevant_proteins_df.index)


        '''
        benchmark_species_lists_directory = "/home/scott/Documents/MATLAB/CAFA2/benchmark/groundtruth/CAFA3/lists"
        benchmark_path = Path(benchmark_species_lists_directory)
        benchmark_species_list_files = benchmark_path.glob(f"{ontology.lower()}_{species_name}*.txt")
        benchmark_species_proteins = set()

        for file in benchmark_species_list_files:
            with open(file) as read_species_proteins_handle:
                proteins = [protein.rstrip() for protein in read_species_proteins_handle.readlines()]
                benchmark_species_proteins.update(proteins)

        benchmark_species_proteins = sorted(list(benchmark_species_proteins))
        '''

        raw_prediction_df = pd.read_csv(
            f,
            engine="python",
            delimiter=prediction_file_delimiter,
            names=("protein", "term", "threshold"),
            skiprows=3,
            skipfooter=1,
            dtype={"protein": "string", "term": "string", "threshold": "float32"},
        )

        # limit the incoming data to the proteins and ontology terms of interest:
        raw_prediction_df = filter_prediction_dataframe(
            raw_prediction_df, proteins=relevant_benchmark_proteins, terms=ontology_terms
        )

        """
        the prediction DataFrame has the following shape with one row
        for every predicted protein/term association: 
        +---------+----------------+------------+-------------+
        |         | protein        | term       |   threshold |
        +=========+================+============+=============+
        | 1534178 | T5592920001831 | GO:0051169 |        0.02 |
        +---------+----------------+------------+-------------+
        | 1534179 | T5592920001831 | GO:0010564 |        0.02 |
        +---------+----------------+------------+-------------+
        | 1534180 | T5592920001831 | GO:0051049 |        0.02 |
        +---------+----------------+------------+-------------+
        | 1534181 | T5592920001831 | GO:0009165 |        0.02 |
        +---------+----------------+------------+-------------+
        | 1534182 | T5592920001831 | GO:0044281 |        0.02 |
        +---------+----------------+------------+-------------+
        """

        raw_prediction_df = get_propagated_dataframe(
            raw_prediction_df, relevant_benchmark_proteins, ontology_dag_df
        )
        """ At this point, the prediction dataframe contains propagated data and has this form:
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        |               |   GO:0000001 |   GO:0000002 |   GO:0000003 |   GO:0000011 |   GO:0000012 |   GO:0000017 |
        +===============+==============+==============+==============+==============+==============+==============+
        | T100900000026 |            0 |            0 |            0 |            0 |            0 |            0 |
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        | T100900000115 |            0 |            0 |            0 |            0 |            0 |            0 |
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        | T100900000116 |            0 |            0 |            0 |            0 |            0 |            0 |
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        | T100900000141 |            0 |            0 |            0 |            0 |            0 |            0 |
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        | T100900000161 |            0 |            0 |            0 |            0 |            0 |            0 |
        +---------------+--------------+--------------+--------------+--------------+--------------+--------------+
        """

        output_filepath = output_path / f"{f.stem}_{ontology_name}.pkl"
        raw_prediction_df.to_pickle(output_filepath)
        print("\t", raw_prediction_df.shape)
        print("\t", datetime.now() - start_t)


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

            '''
            # filter the raw_prediction_df by the benchmark proteins and the benchmark ontology terms:
            # TODO: Refactor the masking here to be less hokey:
            protein_mask = raw_prediction_df['protein'].tolist()
            protein_mask = [v in taxon_ontology_benchmark_df.index for v in protein_mask]
            raw_prediction_df = raw_prediction_df[protein_mask]

            # TODO: This is fragile if the DataFrame columns change in anyway:
            benchmark_terms = taxon_ontology_benchmark_df.columns[2:]
            terms_mask = [v in benchmark_terms for v in raw_prediction_df.loc[:, 'term']]
            raw_prediction_df = raw_prediction_df[terms_mask]
            '''

            '''
            # Create a DataFrame joining the 'raw' prediction data with the DAG propagation DataFrame:
            term_ids = dag_df.columns.values
            raw_prediction_df = pd.merge(raw_prediction_df, dag_df, right_index=True, left_on='term', how='left')
            raw_prediction_df.set_index('protein', inplace=True)

            # identify columns for propagation:
            propagation_mask = raw_prediction_df[term_ids] == 1
            # Here we copy the threshold values to the propagated term columns and
            # simultaneously drop the term and threshold columns. There are still multiple
            # rows per protein:
            raw_prediction_df = raw_prediction_df[term_ids].mask(
                propagation_mask,
                raw_prediction_df["threshold"],
                axis=0
            )

            # Aggregate the rows on protein, keeping the max threshold value for each protein/term pair:
            raw_prediction_df = raw_prediction_df.groupby("protein").aggregate("max")
            '''

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

            # Next, we convert the sparse DataFrame to a dictionary that is more concise
            # and efficient for storage:
            parsed_data = {}

            for protein, term_values in raw_prediction_df.iterrows():
                mask = term_values > 0
                non_zero_annotations = term_values[mask]
                parsed_data[protein] = non_zero_annotations.round(2).to_dict()

            output_filepath = f"./data/ZhangFreddolinoLab/{prediction_file.stem}_{ontology}.json"
            print(f"\tWRITING {output_filepath}")
            with open(output_filepath, 'w') as write_handle:
                json.dump(parsed_data, write_handle, indent=0)

        '''
        parse_raw_predictions(
            predictions_parent_directory,
            benchmark_df_filepath,
            propagation_map_df_filepath,
            ontology_name=ontology,
            output_directory="./data/ZhangFreddolinoLab",
            prediction_file_delimiter="\t",
            species_id_mapping=taxonomy_map
        )
        '''