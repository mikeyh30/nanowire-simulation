import pandas as pd


def update_csv(
    simulation_parameters,
    spectrum_critical_field,
    conductance_critical_field,
    conductance_data_filename,
    spectrum_data_filename,
    conductance_figure_filename,
    spectrum_figure_filename,
    individual_conductance_figure_filename,
    data_file,
):
    df = pd.read_csv(data_file, sep=",")
    simulation_parameters.update(
        {
            "spectrum_critical_field": spectrum_critical_field,
            "conductance_critical_field": conductance_critical_field,
            "conductance_data_filename": conductance_data_filename,
            "spectrum_data_filename": spectrum_data_filename,
            "conductance_figure_filename": conductance_figure_filename,
            "spectrum_figure_filename": spectrum_figure_filename,
            "individual_conductance_figure_filename": individual_conductance_figure_filename,
        },
    )
    newline = df.append(simulation_parameters, ignore_index=True)
    newline.to_csv(data_file, sep=",", index=False)


if __name__ == "__main__":
    df = pd.read_csv("./2020-02-21v12.csv")
    for index, row in df.iterrows():
        update_csv(row.to_dict(), 1, 2, "ee", 4, 5, 6, 7, "./wiresdata.csv")
