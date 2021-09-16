""" This file generates json-formatted versions of CAFA benchmark data from the raw benchmark files.

The json files have this form:

{
    "benchmark_taxon": "CANAX",
    "benchmark_taxon_id": null,
    "benchmark_ontology": "bpo",
    "benchmark_ontology_term_count": 28678,
    "protein_annotations": {
        "T2375610000901": ["GO:0008150", "GO:0030447", "GO:0040007", "GO:0044699"],
        "T2375610001265": ["GO:0008150", "GO:0010570", "GO:0040008", "GO:0045926", ...],
        "T2375610001990": ["GO:0008150", "GO:0030447", "GO:0030448", "GO:0040007", "GO:0044699"],
        ...
    }
}
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from config import taxonomy_map


def main(
    root_benchmark_path_str: str,
    propagation_df_directory_filepath: str,
    #dag_df_directory_filepath: str,
    output_directory_filepath_str: str,
    knowledge_type: int = 1,  # either 1 for partial or 2 for none
    taxonomy_map: dict = {},
    delimiter: str = "\t",
):

    # obo_filepath = "./data/go_cafa3.obo"
    # optional_attrs = ["relationship", "replaced_by", "consider"]
    # optional_relationships = {'part_of', }
    # dag = GODag(obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None)

    root_benchmark_path = Path(root_benchmark_path_str)

    # Benchmark "leaf-only" files are per-ontology, across species
    # each line is a protein/term pair
    benchmark_files = root_benchmark_path.glob("leafonly_*.txt")

    # The "list" files are one file per species/evidence type and contain one protein per line.
    # TODO: Handle types
    species_list_files = list(
        (root_benchmark_path / "lists").glob(f"*_type{knowledge_type}.txt")
    )

    taxonomy_map = {v: k for k, v in taxonomy_map.items()}

    for bf in benchmark_files:
        namespace = bf.stem.split("_")[-1]

        namespace_df = pd.read_csv(bf, delimiter=delimiter, names=("protein", "term"))
        namespace_df.insert(0, "taxon", "")
        namespace_df.insert(1, "taxon_id", "")

        species_list_files_for_namespace = [
            f
            for f in species_list_files
            if namespace.upper() in f.stem.upper() and "_all_" not in f.stem
        ]

        # Loop over the per species/knowledge-type "list" files
        # in order to assign a taxon name and ID to each protein in the raw benchmark:
        for list_file in species_list_files_for_namespace:
            print(f"READING {list_file}")
            namespace, species_short_name, type_str = list_file.stem.split("_")

            with open(list_file, "r") as list_file_handle:
                species_proteins = {line.rstrip() for line in list_file_handle}
                species_mask = namespace_df["protein"].isin(species_proteins)
                namespace_df.loc[species_mask, "taxon"] = species_short_name
                namespace_df.loc[species_mask, "taxon_id"] = taxonomy_map.get(
                    species_short_name
                )

        """
        At this point, we have a DataFrame with this form:
        +------+----------+------------+----------------+------------+
        |      | taxon    | taxon_id   | protein        | term       |
        +====+=========+============+===============+================+
        | 1021 | RAT      | 10116      | T101160003911  | GO:0001540 |
        +------+----------+------------+----------------+------------+
        | 1022 | RAT      | 10116      | T101160003080  | GO:0008134 |
        +------+----------+------------+----------------+------------+
        | 1023 | RAT      | 10116      | T101160003080  | GO:0003682 |
        +------+----------+------------+----------------+------------+
        | 1024 | RAT      | 10116      | T101160003080  | GO:0019213 |
        +------+----------+------------+----------------+------------+
        | 1025 | RAT      | 10116      | T101160005553  | GO:0017124 |
        +------+----------+------------+----------------+------------+
        
        Next, we want to propagate protein->term associations
        """

        propagation_df_filepath = (
            Path(propagation_df_directory_filepath)
            / f"propagation_map_df_{namespace.upper()}.pkl"
        )

        print(f"\tREADING {propagation_df_filepath}")
        propagation_df = pd.read_pickle(propagation_df_filepath)

        taxons = set([taxon for taxon in namespace_df.loc[:, "taxon"] if taxon != ""])

        for taxon in taxons:
            taxon_namespace_df = namespace_df.loc[namespace_df["taxon"] == taxon]
            taxon_namespace_df = taxon_namespace_df.sort_values("protein")

            """ taxon_namespace_df is a species specific subset of the namespace_df
             with this form:
            +------+---------+------------+---------------+------------+
            |      | taxon   |   taxon_id | protein       | term       |
            +======+=========+============+===============+============+
            | 1021 | RAT     |      10116 | T101160003911 | GO:0001540 |
            +------+---------+------------+---------------+------------+
            | 1022 | RAT     |      10116 | T101160003080 | GO:0008134 |
            +------+---------+------------+---------------+------------+
            | 1023 | RAT     |      10116 | T101160003080 | GO:0003682 |
            +------+---------+------------+---------------+------------+
            """

            benchmark_proteins = set(taxon_namespace_df.loc[:, "protein"])
            # Construct a new DataFrame with the species and ontology specific proteins as indices
            # and all the terms in the ontology for columns:
            propagated_benchmark_df = pd.DataFrame(
                index=benchmark_proteins, columns=propagation_df.columns, data=np.nan
            )

            for key, row in taxon_namespace_df.iterrows():
                # get the propagation columns (terms) from the dag DataFrame:
                try:
                    propagation_columns = propagation_df.columns[propagation_df.loc[row.term, :] == 1]
                except KeyError:
                    print(f"COULD NOT PROPAGATE {namespace} TERM {row.term}")

                    continue
                # Populate the DataFrame with the propagation nodes for the given 'leaf' node:
                propagated_benchmark_df.loc[row.protein, propagation_columns] = 1

            # Finally, we are going to merge the simple taxon_namespace_df DataFrame with the propagated_benchmark_df
            taxon_namespace_df = pd.merge(
                taxon_namespace_df,
                propagated_benchmark_df,
                right_index=True,
                left_on="protein",
            )
            taxon_namespace_df.drop("term", axis="columns", inplace=True)
            taxon_namespace_df = taxon_namespace_df.groupby("protein").aggregate("max")
            """ At this point, taxon_namespace_df has this form:
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            | protein       | taxon   |   taxon_id |   GO:0000001 |   GO:0000002 |   GO:0000003 |   GO:0000011 |
            +===============+=========+============+==============+==============+==============+==============+
            | T100900000115 | MOUSE   |      10090 |          nan |          nan |          nan |          ... |
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            | T100900000116 | MOUSE   |      10090 |          nan |          nan |          nan |          ... |
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            | T100900000161 | MOUSE   |      10090 |          nan |          nan |          nan |          ... |
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            | T100900000167 | MOUSE   |      10090 |          nan |          nan |          nan |          ... |
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            | T100900000390 | MOUSE   |      10090 |          nan |          nan |          nan |          ... |
            +---------------+---------+------------+--------------+--------------+--------------+--------------+
            """

            # pickle_filepath = f"./data/benchmark/{namespace.upper()}_{taxon}_{taxonomy_map.get(taxon, '')}_benchmark.pkl"
            # taxon_namespace_df.to_pickle(pickle_filepath)

            benchmark_dict = {
                "benchmark_taxon": taxon,
                "benchmark_taxon_id": taxonomy_map.get(taxon),
                "benchmark_ontology": namespace,
                "benchmark_ontology_term_count": len(propagation_df.index),
                "protein_annotations": {},
            }

            for protein, row in taxon_namespace_df.iterrows():
                row_mask = row.iloc[2:].notna()
                annotated_terms = tuple(row.iloc[2:][row_mask].index)

                benchmark_dict["protein_annotations"][protein] = annotated_terms

            # json_filepath = f"./data/benchmark/{namespace.upper()}_{taxon}_{taxonomy_map.get(taxon, '')}_benchmark.json"
            output_directory_path = Path(output_directory_filepath_str)
            output_directory_path.mkdir(exist_ok=True, parents=True)
            json_filepath = (
                output_directory_path
                / f"{namespace.upper()}_{taxon}_{taxonomy_map.get(taxon, '')}_type_{knowledge_type}_benchmark.json"
            )
            print(f"\tWRITING {json_filepath}")
            with open(json_filepath, "w") as json_write_handle:
                json.dump(benchmark_dict, json_write_handle)


if __name__ == "__main__":
    import sys
    import yaml

    config_filepath = sys.argv[1]

    with open(config_filepath, "r") as config_handle:
        config = yaml.load(config_handle, Loader=yaml.BaseLoader)
        root_benchmark_path_str = config.get("raw_benchmark_path")
        propagation_directory = config.get("propagation_df_directory")

        knowledge_type = config.get("knowledge_type")
        #output_directory_filepath = "./data_v3/parsed_benchmark"
        output_directory_filepath = config.get("benchmark_json_directory")

        main(
            root_benchmark_path_str=root_benchmark_path_str,
            propagation_df_directory_filepath=propagation_directory,
            #dag_df_directory_filepath=dag_df_directory_filepath,
            output_directory_filepath_str=output_directory_filepath,
            knowledge_type=knowledge_type,
            taxonomy_map=taxonomy_map,
            #delimiter=delimiter,
        )
