from pathlib import Path
import pandas as pd

def main():
    ontology = "BPO"
    model_id = 1
    type_id = 1
    # TODO: Address the evaluation type (NK vs PK)
    benchmark_df_filepath = f"../v6/data/benchmarks/benchmark_{ontology}_v2.pkl"
    benchmark_df = pd.read_pickle(benchmark_df_filepath)
    evaluation_directory_filepath = "./data/ZhangFreddolinoLab"
    benchmark_proteins = set(benchmark_df.index)

    benchmark_species_lists_directory = "/home/scott/Documents/MATLAB/CAFA2/benchmark/groundtruth/CAFA3/lists"
    benchmark_path = Path(benchmark_species_lists_directory)
    benchmark_species_list_files = benchmark_path.glob(f"{ontology.lower()}_HUMAN_type{type_id}.txt")

    test = list(benchmark_species_list_files)[0]

    with open(test, 'r') as test_handle:

        test_proteins = set([line.rstrip() for line in test_handle.readlines()])
        benchmark_df = benchmark_df.loc[test_proteins, :]
        print(benchmark_df)

    evaluation_directory_path = Path(evaluation_directory_filepath)
    evaluation_files = evaluation_directory_path.glob(f"*_{model_id}_*_{ontology}*")

    for ef in evaluation_files:

        if '9606' not in str(ef):
            continue

        prediction_df = pd.read_pickle(ef)
        print(set(prediction_df.index))
        print(set(test_proteins))
        assert set(prediction_df.index) == set(test_proteins)

        #prediction_df = prediction_df.loc[test_proteins, :]
        #print("++++++++++++++++++++++++++++++")
        #print(prediction_df)


if __name__ == "__main__":
    main()
