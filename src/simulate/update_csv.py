import pandas as pd
import h5py

#def update_csv(
#    simulation_parameters,
#    data_file,
#):
#    df = pd.read_csv(data_file, sep=",")
#    newline = df.append(simulation_parameters, ignore_index=True)
#    newline.to_csv(data_file, sep=",", index=False)
#==========
# def update_csv(
#     simulation_parameters,
#     spectrum_critical_field,
#     conductance_critical_field,
#     conductance_data_filename,
#     spectrum_data_filename,
#     conductance_figure_filename,
#     spectrum_figure_filename,
#     individual_conductance_figure_filename,
#     data_file,
# ):
#     df = pd.read_csv(data_file, sep=",")
#     simulation_parameters.update(
#         {
#             "spectrum_critical_field": spectrum_critical_field,
#             "conductance_critical_field": conductance_critical_field,
#             "conductance_data_filename": conductance_data_filename,
#             "spectrum_data_filename": spectrum_data_filename,
#             "conductance_figure_filename": conductance_figure_filename,
#             "spectrum_figure_filename": spectrum_figure_filename,
#             "individual_conductance_figure_filename": individual_conductance_figure_filename,
#         },
#     )
#     newline = df.append(simulation_parameters, ignore_index=True)
#     newline.to_csv(data_file, sep=",", index=False)

def update_csv(iteration, spectrum_critical_field, conductance_critical_field, data_file):
    df = pd.read_csv(data_file, sep=",")
    df.at[iteration, 'spectrum_critical_field'] = spectrum_critical_field
    df.at[iteration, 'conductance_critical_field'] = conductance_critical_field
    df.to_csv(data_file, sep=",", index=False, mode="w")


def add_dataset_hdf(groupname, data_file, **kwargs):
    with h5py.File(data_file,'a') as file:
        grp=file[groupname]
        for key, value in kwargs.items():
            grp.create_dataset(key,data=value)


import time

if __name__ == "__main__":
    df = pd.read_csv("/home/michael/Documents/UCL/hydrogen/scratch/testy2/testy2.csv")
    start = time.time()
    for index, row in df.iterrows():
        update_csv(index, index, 1, "/home/michael/Documents/UCL/hydrogen/scratch/testy2/testy2.csv")
    end = time.time()
    print(end-start)
