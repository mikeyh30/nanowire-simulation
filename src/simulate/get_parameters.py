from numpy import arange
import os
import yaml


def get_yml(yml_file_name):
    yml_file_name = os.path.join(os.path.dirname(__file__), "../../", yml_file_name)
    with open(yml_file_name) as yml_file:
        yml_file_entries = yaml.load(yml_file, Loader=yaml.FullLoader)
    return yml_file_entries


def get_scratch():
    scratch = get_yml("globals.yml")["directories"]["scratch"]
    return scratch


def get_simulation_parameters():
    use_units = get_yml("globals.yml")["use_units"]
    if use_units:
        return get_yml("sim_parameters_units.yml")
    else:
        return get_yml("sim_parameters_unitless.yml")


if __name__ == "__main__":
    print(get_simulation_parameters())
