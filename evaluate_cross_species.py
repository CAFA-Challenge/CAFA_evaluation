from pathlib import Path
import numpy as np
import pandas as pd


if __name__ == "__main__":


    root_path = Path("./data/working/ZhangFreddolinoLab/")
    species_specific_files = list(root_path.glob("*.pkl"))

    # Collect all thresholds across all of the species predictions:
    all_thresholds = set()

    for file in species_specific_files:
        df = pd.read_pickle(file)
        thresholds = df.index.get_level_values(1)

        all_thresholds.update(thresholds)




    all_thresholds = sorted(list(all_thresholds))
    columns = ('threshold', 'ontology', 'taxons', 'average_precision', 'average_recall', 'average_weighted_precision', 'average_weighted_recall')

    # Construct a dataframe to hold the computed results:
    matrix = np.zeros((len(all_thresholds), len(columns)))
    results_df = pd.DataFrame(
        data=matrix,
        columns=columns
    )
    results_df = results_df.astype({'ontology': 'str'})
    results_df['threshold'] = all_thresholds
    results_df.set_index('threshold', drop=True, inplace=True)


    for threshold in all_thresholds:
        all_species_at_threshold_df = None

        for file in species_specific_files:

            df = pd.read_pickle(file)

            threshold_mask = df.index.get_level_values(1) >= threshold
            df = df[threshold_mask]
            ''' 
            try:
                df = df.loc[pd.IndexSlice[:, threshold], :]
            except KeyError:
                # the threshold at hand does not exist in the species at hand
                continue
            '''
            if all_species_at_threshold_df is None:
                all_species_at_threshold_df = df
            else:
                all_species_at_threshold_df = all_species_at_threshold_df.append(df)


        protein_count = all_species_at_threshold_df.shape[0]

        average_precision = all_species_at_threshold_df.loc[:, 'precision'].sum() / protein_count
        average_recall = all_species_at_threshold_df.loc[:, 'recall'].sum() / protein_count
        average_weighted_precision = all_species_at_threshold_df.loc[:, 'weighted_precision'].sum() / protein_count
        average_weighted_recall = all_species_at_threshold_df.loc[:, 'weighted_recall'].sum() / protein_count

        results_df.loc[threshold, 'ontology'] = all_species_at_threshold_df.iloc[0].ontology
        results_df.loc[threshold, 'average_precision'] = average_precision
        results_df.loc[threshold, 'average_recall'] = average_recall
        results_df.loc[threshold, 'average_weighted_precision'] = average_weighted_precision
        results_df.loc[threshold, 'average_weighted_recall'] = average_weighted_recall

        taxons = ",".join(set(all_species_at_threshold_df.loc[:, 'taxon']))
        results_df.loc[threshold, 'taxons'] = taxons


    print(results_df.to_markdown(tablefmt='grid'))






    import sys
    sys.exit(0)



    filepath = "./data/working/ZhangFreddolinoLab/HUMAN_cco_1.pkl"
    df = pd.read_pickle(filepath)

    for row in df.iterrows():
        print(row)


    all_thresholds = set(df.index.get_level_values(1))
    all_thresholds = sorted(list(all_thresholds))

    for threshold in all_thresholds:
        print("\n===========================")
        print(threshold)
        print(df.loc[pd.IndexSlice[:, threshold], :])
