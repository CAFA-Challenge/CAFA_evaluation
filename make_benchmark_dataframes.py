from pathlib import Path
import numpy as np
import pandas as pd

delimiter = "\t"
root_benchmark_path_str = "/home/scott/Documents/MATLAB/CAFA2/benchmark/groundtruth/CAFA3/"
root_benchmark_path = Path(root_benchmark_path_str)

# Benchmark "leaf-only" files are per-ontology, across species
# each line is a protein/term pair
benchmark_files = root_benchmark_path.glob("leafonly_*.txt")

# The "list" files are one file per species/evidence type and contain one protein per line.
species_list_files = list((root_benchmark_path / "lists").glob("*_type1.txt"))

# The 'list' files use non-standard short names, while the submitted predictions use official
# taxon IDs, so we need a way to translate between the two:
species_map = {
    9606: 'HUMAN',
    3702: 'ARATH',  # Arabidopsis
    7227: 'DROME',  # Drosophila melanogaster
    10090: 'MOUSE',  # Mus musculus
    10116: 'RAT',  # Rattus norvegicus
}

species_map = {v:k for k, v in species_map.items()}

for bf in benchmark_files:
    namespace = bf.stem.split("_")[-1]

    namespace_df = pd.read_csv(bf, delimiter=delimiter, names=('protein', 'term'))
    namespace_df.insert(0, 'taxon', '')
    namespace_df.insert(1, 'taxon_id', '')

    species_list_files_for_namespace = [f for f in species_list_files if namespace.upper() in f.stem.upper() and '_all_' not in f.stem]

    # Loop over the per species/knowledge-type "list" files
    # in order to assign a taxon name and ID to each protein in the raw benchmark:
    for list_file in species_list_files_for_namespace:
        namespace, species_short_name, type_str = list_file.stem.split("_")

        with open(list_file, 'r') as list_file_handle:
            species_proteins = {line.rstrip() for line in list_file_handle}
            species_mask = namespace_df['protein'].isin(species_proteins)
            namespace_df.loc[species_mask, 'taxon'] = species_short_name
            namespace_df.loc[species_mask, 'taxon_id'] = species_map.get(species_short_name, '')

    '''
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
    '''

    dag_df_filepath = f"../v6/data/propagation/propagation_map_df_{namespace.upper()}.pkl"
    print(f"\tREADING {dag_df_filepath}")
    dag_df = pd.read_pickle(dag_df_filepath)

    taxons = set([taxon for taxon in namespace_df.loc[:, 'taxon'] if taxon != ''])

    for taxon in taxons:
        taxon_namespace_df = namespace_df.loc[namespace_df['taxon'] == taxon]
        taxon_namespace_df = taxon_namespace_df.sort_values('protein')

        ''' taxon_namespace_df is a species specific subset of the namespace_df
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
        '''

        benchmark_proteins = set(taxon_namespace_df.loc[:, 'protein'])
        # Construct a new DataFrame with the species and ontology specific proteins as indices
        # and all the terms in the ontology for columns:
        propagated_benchmark_df = pd.DataFrame(index=benchmark_proteins, columns=dag_df.columns, data=np.nan)

        for key, row in taxon_namespace_df.iterrows():
            # get the propagation columns (terms) from the dag DataFrame:
            try:
                propagation_columns = dag_df.columns[dag_df.loc[row.term, :] == 1]
            except KeyError:
                print(f"COULD NOT PROPAGATE {namespace} TERM {row.term}")
                continue
            # Populate the DataFrame with the propagation nodes for the given 'leaf' node:
            propagated_benchmark_df.loc[row.protein, propagation_columns] = 1

        # Finally, we are going to merge the simple taxon_namespace_df DataFrame with the propagated_benchmark_df
        taxon_namespace_df = pd.merge(taxon_namespace_df, propagated_benchmark_df, right_index=True, left_on='protein') #left_on='protein', right_on='index')
        taxon_namespace_df.drop('term', axis='columns', inplace=True)
        taxon_namespace_df = taxon_namespace_df.groupby("protein").aggregate("max")
        ''' At this point, taxon_namespace_df has this form:
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
        '''

        pickle_filepath = f"./data/benchmark/{namespace.upper()}_{taxon}_{species_map.get(taxon, '')}_benchmark.pkl"
        taxon_namespace_df.to_pickle(pickle_filepath)
