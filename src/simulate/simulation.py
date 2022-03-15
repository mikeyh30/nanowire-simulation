import matplotlib.pyplot as plt
import numpy as np
import pickle
from nanowire.nanowire import Nanowire
from simulate.update_csv import update_csv
from simulate.get_parameters import get_scratch, get_yml
import argparse
import pandas as pd


def save_model_figure(nanowire, suffix, data_folder):
    fig = plt.figure()
    fig.suptitle("Model grid, and nanomagnet fields")
    ax_model = fig.add_subplot(3, 1, 1)
    ax_model.axis("equal")
    ax_x = fig.add_subplot(3, 1, 2)
    ax_y = fig.add_subplot(3, 1, 3)
    nanowire.plot(ax_model, ax_x, ax_y)
    ax_model.set_ylabel(r"Length ($\AA$)")
    ax_x.set_ylabel(r"$M_x$")
    ax_y.set_ylabel(r"$M_y$")
    ax_y.set_xlabel(r"Length ($\AA$)")
    xlim = fig.gca().figure.axes[0].dataLim.intervalx
    plt.xlim(xlim)
    plt.savefig(data_folder + "/modelfig/" + suffix + ".eps", format="eps", dpi=1000)
    plt.close()


def spectrum(spectrum_data, suffix, data_folder):
    pickle.dump(spectrum_data, open(data_folder + "/spec/spec_" + suffix + ".dat", "wb"))

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (7, 5)
    ax = fig.gca()
    ax.plot(spectrum_data["B"], spectrum_data["E"])
    ax.set_xlabel("Zeeman Field Strength (T)")
    ax.set_ylabel("Energies (eV)")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    fig.savefig(data_folder + "/fig-spectrum/model" + suffix + ".png")
    plt.close(fig)
    return spectrum_data["CritB"]


def density(density_data, suffix, data_folder, params):
    pickle.dump(density_data, open(data_folder + "/dens/dens_" + suffix + ".dat", "wb"))
    wf1=[]
    wf2=[]
    for x in range(params['wire_width']*params['wire_length']):
        wf1.append(density_data[x][0])
        wf2.append(density_data[x][1])
    wf = np.array(wf1)+np.array(wf2)
    wfplot=wf.reshape(params['wire_length'],params['wire_width']).transpose()

    mult = params['hopping_distance']/10 # divide by ten for A->nm
    xlist = np.linspace(0, params['wire_length']*mult, params['wire_length'])
    ylist = np.linspace(0, params['wire_width']*mult, params['wire_width'])
    x, y = np.meshgrid(xlist, ylist)  
    fig,ax=plt.subplots(1,1,figsize=(15,4))

    cp = ax.contourf(x, y, wfplot,50)
    fig.colorbar(cp,orientation="horizontal", label='Density') # Add a colorbar to a plot
    ax.axis('scaled')
    ax.set_ylabel('Size (nm)', fontsize=12)
    fig.savefig(data_folder + "/fig-dens/model" + suffix + ".eps", format="eps", dpi=1000)
    fig.close()


def magnetization_spectrum(spectrum_data, suffix, data_folder):
    pickle.dump(spectrum_data, open(data_folder + "/spec/spec_" + suffix + ".dat", "wb"))

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (7, 5)
    ax = fig.gca()
    ax.plot(spectrum_data["M"], spectrum_data["E"])
    ax.set_xlabel("Magnetization (T)")
    ax.set_ylabel("Energies (eV)")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    fig.savefig(data_folder + "/fig-spectrum/model" + suffix + ".png")
    plt.close(fig)
    return spectrum_data["CritM"]


def conductance(conductance_data, suffix, data_folder):
    pickle.dump(conductance_data, open(data_folder + "/cond/cond_" + suffix + ".dat", "wb"))

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (8, 5)
    ax = fig.gca()
    contour = ax.contourf(
        conductance_data["B"],
        conductance_data["BiasV"],
        conductance_data["Cond"],
        100,
        cmap="viridis",
    )
    ax.set_xlabel("Zeeman Field Strength (T)")
    ax.set_ylabel("Bias (V)")
    cbar = plt.colorbar(contour)
    cbar.ax.set_ylabel("Conductance [e^2/h]")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    fig.savefig(data_folder + "/fig-conductance/model" + suffix + ".png")
    plt.close(fig)
    return conductance_data["CritB"]


def individual_conductance(data, suffix, data_folder, index_slice=30):
    plt.rcParams["figure.figsize"] = (7, 5)

    cond = np.transpose(data["Cond"])
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(data["BiasV"], cond[index_slice])
    ax.set_xlabel("Bias (V)")
    ax.set_ylabel("Conductance [e^2/h]")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    fig.savefig(data_folder + "/fig-ind-conductance/model" + suffix + ".png")
    plt.close(fig)


def simulation_single(
    params,
    row,
    date="no-date",
    scratch="./Scratch/",
    simulate_conductance=True,
    simulate_spectrum=True,
    simulate_magnetization_spectrum=True,
    simulate_density=False,
):

    data_suffix = "simulation{}".format(row)

    params.update({"row": row})

    nanowire = Nanowire(params)

    data_folder = scratch + date

    save_model_figure(nanowire, data_suffix, data_folder)

    if simulate_spectrum:
        # Generate spectrum data and figure
        spectrum_data = nanowire.spectrum(B_values=np.linspace(0, params["b_max"], 81))
        spectrum_critical_field = spectrum(spectrum_data, data_suffix, data_folder)
        params.update({"spectrum_critical_field": spectrum_critical_field})
        params.update({"superconducting_gap": spectrum_data["superconducting_gap"]})

    if simulate_magnetization_spectrum:
        # Generate spectrum data and figure
        spectrum_data = nanowire.magnetization_spectrum(M_values=np.linspace(0, params["m_max"], 81))
        mag_spectrum_critical_field = magnetization_spectrum(spectrum_data, data_suffix, data_folder)
        params.update({"mag_spectrum_critical_field": mag_spectrum_critical_field})

    if simulate_conductance:
        # Generate data of spectrum and conductance. This takes time
        t = nanowire.parameters["t"]
        energies = np.arange(-0.120 * t, 0.120 * t, 0.001 * t)
        # Generate conductance data and figure
        conductance_data = nanowire.conductances(B_values=np.linspace(0, params["b_max"], 81), energies=energies)
        conductance_critical_field = conductance(conductance_data, data_suffix, data_folder)

        # Save figure of the conductance at a given field.
        individual_conductance(conductance_data, data_suffix, data_folder)
        params.update({"conductance_critical_field": conductance_critical_field})

    if simulate_density:
        density_data = nanowire.density(2,4)
        density(density_data,data_suffix,data_folder,params)

    # Log which data has been saved.
    update_csv(
        params,
        data_folder + "/wiresdata.csv",
    )


def simulation_all(params, row, date="no-date", scratch="./Scratch/"):
    new_params = params
    for no_magnets in params["Ns"]:
        new_params["no_magnets"] = no_magnets
        simulation_single(new_params, row, date, scratch)


def simulation_all_csv(
    csv_file, date, scratch, simulate_conductance, simulate_spectrum, simulate_magnetization_spectrum,simulate_density
):
    df = pd.read_csv(csv_file)
    for index, row in df.iterrows():
        simulation_single(
            row.to_dict(),
            row=f"{index:0>5}",
            date=date,
            scratch=scratch,
            simulate_conductance=simulate_conductance,
            simulate_magnetization_spectrum=simulate_magnetization_spectrum,
            simulate_spectrum=simulate_spectrum,
            simulate_density=simulate_density,
        )


def main():
    parser = argparse.ArgumentParser(description="take the csv, and the line number")
    parser.add_argument("csv_file", metavar="filename", type=str)
    parser.add_argument("date", type=str)

    args = parser.parse_args()

    simulate_conductance = get_yml("globals.yml")["simulate_conductance"]
    simulate_magnetization_spectrum = get_yml("globals.yml")["simulate_magnetization_spectrum"]
    simulate_spectrum = get_yml("globals.yml")["simulate_spectrum"]
    simulate_density = get_yml("globals.yml")["simulate_density"]

    simulation_all_csv(
        args.csv_file,
        args.date,
        get_scratch(),
        simulate_conductance,
        simulate_spectrum,
        simulate_magnetization_spectrum,
        simulate_density,
    )


if __name__ == "__main__":
    main()
