import pandas as pd


def update_csv(simulation_parameters, data_file):
    df = pd.read_csv(data_file, sep=",")
    # newline = df.append(simulation_parameters, ignore_index=True)
    newline2 = pd.concat([df, pd.DataFrame(simulation_parameters, index=[0])])
    # newline.to_csv(data_file, sep=",", index=False)
    newline2.to_csv(data_file, sep=",", index=False)


if __name__ == "__main__":
    pass
