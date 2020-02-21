from numpy import arange
import os
import yaml
from simulate.get_parameters import get_simulation_parameters

# from simulate.setup_sim import 

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
):
    df = pd.read_csv("./wiresdata.csv", sep=",")
    print("uu",type(simulation_parameters))
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
    newline = df.append(simulation_parameters,ignore_index=True)
    newline.to_csv("./wiresdata.csv", sep=",", index=False)

if __name__ == "__main__":
    df = pd.read_csv("./2020-02-21v12.csv")
    for index, row in df.iterrows():
        print("ee",type(row))
        print("oo",type(get_simulation_parameters()))
        update_csv(row.to_dict(),1,2,3,4,5,6,7)


