import pandas as pd


if __name__ == "__main__":

    filepath = "./data/working/ZhangFreddolinoLab/HUMAN_cco_1.pkl"

    df = pd.read_pickle(filepath)

    for row in df.iterrows():
        print(row)


    threshold_set = set(df.index.get_level_values(1))
    print(threshold_set)

    for threshold in threshold_set:
        print("\n===========================")
        print(threshold)
        print(df.loc[pd.IndexSlice[:, threshold], :])
