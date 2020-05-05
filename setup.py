import os
from shutil import copyfile
from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop


def run_installer(command_subclass):
    command_subclass_run = command_subclass.run

    def modified_run(self):
        # Complete traditional installation
        command_subclass_run(self)
        # Define config file locations
        globals_sample = os.path.join(os.getcwd(), "config", "globals-sample.yml")
        units_sample = os.path.join(os.getcwd(), "config", "sim_parameters_units-sample.yml")
        unitless_sample = os.path.join(os.getcwd(), "config", "sim_parameters_unitless-sample.yml")
        globals_file = os.path.join(os.getcwd(), "globals.yml")
        units_file = os.path.join(os.getcwd(), "sim_parameters_units.yml")
        unitless_file = os.path.join(os.getcwd(), "sim_parameters_unitless.yml")
        # Copy config files across
        copyfile(globals_sample, globals_file)
        copyfile(units_sample, units_file)
        copyfile(unitless_sample, unitless_file)

        # Update globals
        # Read in the file
        with open('globals.yml', 'r') as file:
            filedata = file.read()

        scratch = os.path.join(os.getcwd(), "data/")
        src = os.path.join(os.getcwd(), "src/")
        hpc_control = os.path.join(os.getcwd(), "hpc-control/")

        # Replace the target string
        filedata = filedata.replace("/path/to/scratch/folder/", scratch)
        filedata = filedata.replace("/path/to/src/folder/", src)
        filedata = filedata.replace("/path/to/hpc-control/folder/", hpc_control)

        # Write the file out again
        with open('globals.yml', 'w') as file:
            file.write(filedata)

    command_subclass.run = modified_run
    return command_subclass


@run_installer
class SetupConfigurationInstall(install):
    pass

@run_installer
class SetupConfigurationDevelop(develop):
    pass


setup(
    name="nanowire-simulation",
    author="Michael Hynes",
    version="0.2.0",
    packages=find_packages("src", exclude=[]),
    python_requires=">=3.5",
    package_dir={"": "src"},
    install_requires=[
        "kwant>=1.4.1",
        "numpy>=1.11.0",
        "scipy>=0.17.0",
        "matplotlib>=1.5.1",
        "sympy>=0.7.6",
        "tinyarray>=1.2",
        "tqdm",
        "pyyaml",
        "pandas",
        "pre-commit",
    ],
    entry_points={
        "console_scripts": [
            "nwsetup = simulate.setup_sim:main",
            "simulate = simulate.simulation:main",
            "hpc_simulate = simulate.hpc_simulation:main",
        ]
    },
    cmdclass={
        "install": SetupConfigurationInstall,
        "develop": SetupConfigurationDevelop,
    },
)
