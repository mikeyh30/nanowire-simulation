import kwant
import matplotlib.pyplot as plt
import numpy as np
import pickle
from nanowire import Nanowire
from update_csv import update_csv
from emailer import send_finished_script_email
from simulation_parameters import simulation_parameters
import argparse
import os

def save_model_figure(nanowire, suffix, data_folder):
    ax = nanowire.plot()
    ax.savefig(data_folder + "/modelfig/" + suffix+ ".png")
    plt.close()

def spectrum(spectrum_data, suffix, data_folder):
    pickle.dump(spectrum_data, open(data_folder + "/spec/spec_" + suffix + ".dat", "wb"))
    
    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (7,5)
    ax = fig.gca()
    ax.plot(spectrum_data["B"], spectrum_data["E"])
    ax.set_xlabel("Zeeman Field Strength [B]")
    ax.set_ylabel("Energies [t]")
    fig.savefig(data_folder + "/fig-spectrum/model" + suffix + ".png")

    return spectrum_data["CritB"]

def conductance(conductance_data, suffix, data_folder):
    pickle.dump(conductance_data, open(data_folder + "/cond/cond_" + suffix + ".dat", "wb"))

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (8,5)
    ax = fig.gca()
    contour = ax.contourf(conductance_data["B"], conductance_data["BiasV"], conductance_data["Cond"], 100, cmap="viridis")
    ax.set_xlabel("Zeeman Field Strength [B]")
    ax.set_ylabel("Bias V [t]")
    cbar = plt.colorbar(contour)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    fig.savefig(data_folder + "/fig-conductance/model" + suffix + ".png")

    return conductance_data["CritB"]

def individual_conductance(data, suffix, data_folder, index_slice=30):
    plt.rcParams["figure.figsize"] = (7,5)
    
    cond = np.transpose(data["Cond"])
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(data["BiasV"], cond[index_slice])
    ax.set_xlabel("Bias V [t]")
    ax.set_ylabel("Conductance [e^2/h]")
    fig.savefig(data_folder + "/fig-ind-conductance/model" + suffix + ".png")
    plt.close()

def simulation_single(params,row='skip'):
    if row == 'skip':
        data_suffix = "w{0}_no{1}_eM{2:3.2f}_mu{3}_al{4}_M{5:4.2f}_added{6}_ratio{7:4.2f}".format( 
        params['wire_width'], params['N'], params['effective_mass'], params['muSc'],
        params['alpha'], params['M'], int(params['added_sinusoid']), params['ratio'])
    else:
        data_suffix = "simulation{}".format(row)

    nanowire = Nanowire(width=params['wire_width'],
                        noMagnets=params['N'],
                        effective_mass=params['effective_mass'],
                        muSc=params['muSc'],
                        alpha=params['alpha'],
                        M=params['M'],
                        addedSinu=params['added_sinusoid'],
                        stagger_ratio=params['ratio'],
                        mu=params['mu'],
                        delta=params['delta'],
                        barrier=params['barrier'])

    data_folder = "/home/ucapmhy/Scratch/2019-11-15"
    os.makedirs(data_folder + '/modelfig',exist_ok=True)
    os.makedirs(data_folder + '/cond',exist_ok=True)
    os.makedirs(data_folder + '/spec',exist_ok=True)
    os.makedirs(data_folder + '/fig-conductance',exist_ok=True)
    os.makedirs(data_folder + '/fig-ind-conductance',exist_ok=True)
    os.makedirs(data_folder + '/fig-spectrum',exist_ok=True)
    
    save_model_figure(nanowire, data_suffix, data_folder)

    # Generate data of spectrum and conductance. This takes time
    energies = np.arange(-0.120*nanowire.t, 0.120*nanowire.t, 0.001*nanowire.t)

    spectrum_data = nanowire.spectrum(bValues=np.linspace(0, params['b_max'], 81))
    conductance_data = nanowire.conductances(bValues=np.linspace(0, params['b_max'], 81),
                                             energies=energies)

    # Save figures and data, and get critical fields.
    spectrum_critical_field = spectrum(spectrum_data, data_suffix, data_folder)
    conductance_critical_field = conductance(conductance_data, data_suffix, data_folder)

    # Save figure of the conductance at a given field.
    individual_conductance(conductance_data, data_suffix, data_folder)

    # Filenames of the saved data and figures.
    conductance_data_filename = data_folder + "/cond/cond_" + data_suffix + ".dat"
    spectrum_data_filename = data_folder + "/spec/spec_" + data_suffix + ".dat"
    conductance_figure_filename = data_folder + "/fig-conductance/model" + data_suffix + ".png"
    spectrum_figure_filename = data_folder + "/fig-spectrum/model" + data_suffix + ".png"
    individual_conductance_figure_filename = data_folder + "/fig-ind-conductance/model" + data_suffix + ".png"

    # Log which data has been saved. Alter this function.
    update_csv(params['wire_width'],
               params['N'],
               params['effective_mass'],
               params['muSc'],
               params['alpha'],
               params['M'],
               params['added_sinusoid'],
               params['ratio'],
               conductance_data_filename,
               spectrum_data_filename,
               conductance_figure_filename,
               spectrum_figure_filename,
               individual_conductance_figure_filename,
               '/home/ucapmhy/Scratch/2019-11-15/wiresdata.csv',
               spectrum_critical_field,
               conductance_critical_field,
               data_folder
               )

def simulation_all(params):
    new_params = params
    for N in params['Ns']:
        new_params['N'] = N
        simulation_single(new_params)


if __name__ == "__main__":
    simulation_all(simulation_parameters)
