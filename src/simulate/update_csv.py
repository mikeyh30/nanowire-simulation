import pandas as pd


def update_csv(
    width,
    noMagnets,
    effectiveMass,
    muSc,
    alpha_R,
    M,
    added_sinusoid,
    ratio,
    conductance_data_filename,
    spectrum_data_filename,
    conductance_figure_filename,
    spectrum_figure_filename,
    individual_conductance_figure_filename,
    csv_filename,
    spectrum_critical_field,
    conductance_critical_field,
    data_folder,
):
    df = pd.read_csv(data_folder + "/wiresdata.csv", sep=",")

    newline = df.append(
        {
            "wire_width": width,
            "no_magnets": noMagnets,
            "effective_mass": effectiveMass,
            "muSc": muSc,
            "alpha_R": alpha_R,
            "M": M,
            "added_sinusoid": added_sinusoid,
            "ratio": ratio,
            "conductance_data_filename": conductance_data_filename,
            "spectrum_data_filename": spectrum_data_filename,
            "conductance_figure_filename": conductance_figure_filename,
            "spectrum_figure_filename": spectrum_figure_filename,
            "individual_conductance_figure_filename": individual_conductance_figure_filename,
            "spectrum_critical_field": spectrum_critical_field,
            "conductance_critical_field": conductance_critical_field,
        },
        ignore_index=True,
    )

    newline.to_csv(data_folder + "/wiresdata.csv", sep=",", index=False)


if __name__ == "__main__":
    update_csv(
        1, 2, 3, 4, 5, 6, True, 123, 12, 34, 56, 32, 24, "data/wiresdata.csv", 1, 1
    )
