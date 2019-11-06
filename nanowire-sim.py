import kwant
import matplotlib.pyplot as plt
import numpy as np
import pickle
from nwObjects import Nanowire
from update_csv import update_csv
from emailer import send_finished_script_email

#------------Simulation parameters-----------------
minN = 6
maxN = 9

simulation_parameters = dict(
    wire_width = 7,
    Ns = np.arange(minN,maxN+1,1), #N is the number of magnets
    ratio = 0.5, #ratios=[i for i in np.arange(0.20,0.50,0.05)], # relative ratios of nanomagnet widths.
    M = 0.1, #B field strength from nanomagnets.
    added_sinusoid = True, # Indicates presence of nanomagnets
    ## SOI terms ##
    effective_mass=.5, #Of electrons
    muSc=.22, #Chemical potential in the nanowire.
    alpha=.0 #Rashba parameter
)
#--------------------------------------------------

def save_model_figure(nanowire, suffix):
    ax = nanowire.plot()
    ax.savefig("data/modelfig/model" + suffix+ ".png")
    plt.close()

def spectrum(spectrum_data, suffix):
    pickle.dump(spectrum_data, open("data/spec/spec_" + suffix + ".dat", "wb"))
    
    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (7,5)
    ax = fig.gca()
    ax.plot(spectrum_data["B"], spectrum_data["E"])
    ax.set_xlabel("Zeeman Field Strength [B]")
    ax.set_ylabel("Energies [t]")
    fig.savefig("data/fig-spectrum/model" + suffix + ".png")

    return spectrum_data["CritB"]

def conductance(conductance_data, suffix):
    pickle.dump(conductance_data, open("data/cond/cond_" + suffix + ".dat", "wb"))

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (8,5)
    ax = fig.gca()
    contour = ax.contourf(conductance_data["B"], conductance_data["BiasV"], conductance_data["Cond"], 100, cmap="viridis")
    ax.set_xlabel("Zeeman Field Strength [B]")
    ax.set_ylabel("Bias V [t]")
    cbar = plt.colorbar(contour)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    fig.savefig("data/fig-conductance/model" + suffix + ".png")

    return conductance_data["CritB"]

def individual_conductance(data, suffix, index_slice=30):
    plt.rcParams["figure.figsize"] = (7,5)
    
    cond = np.transpose(data["Cond"])
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(data["BiasV"], cond[index_slice])
    ax.set_xlabel("Bias V [t]")
    ax.set_ylabel("Conductance [e^2/h]")
    fig.savefig("data/fig-ind-conductance/model" + suffix + ".png")
    plt.close()

def simulation_single(params):
    data_suffix = "w{0}_no{1}_eM{2:3.2f}_mu{3}_al{4}_M{5:4.2f}_added{6}_ratio{7:4.2f}".format( 
        params['wire_width'], params['N'], params['effective_mass'], params['muSc'],
        params['alpha'], params['M'], int(params['added_sinusoid']), params['ratio'])

    nanowire = Nanowire(width=params['wire_width'],
                        noMagnets=params['N'],
                        effective_mass=params['effective_mass'],
                        muSc=params['muSc'],
                        alpha=params['alpha'],
                        M=params['M'],
                        addedSinu=params['added_sinusoid'],
                        stagger_ratio=params['ratio'])

    save_model_figure(nanowire, data_suffix)

    # Generate data of spectrum and conductance. This takes time
    spectrum_data = nanowire.spectrum(bValues=np.linspace(0, .4, 81))
    conductance_data = nanowire.conductances(bValues=np.linspace(0, .4, 81))

    # Save figures and data, and get critical fields.
    spectrum_critical_field = spectrum(spectrum_data, data_suffix)
    conductance_critical_field = conductance(conductance_data, data_suffix)

    # Save figure of the conductance at a given field.
    individual_conductance(conductance_data, data_suffix)

    # Filenames of the saved data and figures.
    conductance_data_filename = "data/cond/cond_" + data_suffix + ".dat"
    spectrum_data_filename = "data/spec/spec_" + data_suffix + ".dat"
    conductance_figure_filename = "data/fig-conductance/model" + data_suffix + ".png"
    spectrum_figure_filename = "data/fig-spectrum/model" + data_suffix + ".png"
    individual_conductance_figure_filename = "data/fig-ind-conductance/model" + data_suffix + ".png"

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
               'data/wiresdata.csv',
               spectrum_critical_field,
               conductance_critical_field
               )

def simulation_all(params):
    new_params = params
    for N in params['Ns']:
        new_params['N'] = N
        simulation_single(new_params)
    send_finished_script_email()


if __name__ == "__main__":
    simulation_all(simulation_parameters)
