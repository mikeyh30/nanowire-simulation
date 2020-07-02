import pandas as pd


def update_csv(simulation_parameters, data_file):
    df = pd.read_csv(data_file, sep=",")
    newline = df.append(simulation_parameters, ignore_index=True)
    newline.to_csv(data_file, sep=",", index=False)


if __name__ == "__main__":
    df = pd.read_csv("./2020-02-21v12.csv")
    for index, row in df.iterrows():
        update_csv(row.to_dict(), 1, 2, "ee", 4, 5, 6, 7, "./wiresdata.csv")
