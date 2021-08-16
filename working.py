from pathlib import Path
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


    for threshold in all_thresholds:
        all_species_at_threshold_df = None

        for file in species_specific_files:

            df = pd.read_pickle(file)

            try:
                df = df.loc[pd.IndexSlice[:, threshold], :]
            except KeyError:
                # the threshold at hand does not exist in the species at hand
                continue

            if all_species_at_threshold_df is None:
                all_species_at_threshold_df = df
            else:
                all_species_at_threshold_df = all_species_at_threshold_df.append(df)


        print(threshold)
        print(all_species_at_threshold_df.loc[:, ('ontology', 'taxon_id', 'taxon', 'weighted_precision', 'weighted_recall', 'precision', 'recall', 'tp', 'tn', 'fn', 'fp')].to_markdown(tablefmt='grid'))
        protein_count = all_species_at_threshold_df.shape[0]

        average_precision = all_species_at_threshold_df.loc[:, 'precision'].sum() / protein_count
        average_recall = all_species_at_threshold_df.loc[:, 'recall'].sum() / protein_count

        print(f"AVERAGE PRECISION: {average_precision}, AVERAGE RECALL: {average_recall}")
        average_weighted_precision = all_species_at_threshold_df.loc[:, 'weighted_precision'].sum() / protein_count
        average_weighted_recall = all_species_at_threshold_df.loc[:, 'weighted_recall'].sum() / protein_count
        print(f"AVERAGE WEIGHTED PRECISION: {average_weighted_precision}, AVERAGE WEIGHTED RECALL: {average_weighted_recall}")


        print("")




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
