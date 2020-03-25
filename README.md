# Simulations of transport properties of a nanowire in a magnetic field

This repository contains code for performing simulations of an InAs nanowire using the kwant package. It is a continuation of the work done in <https://github.com/dk1713/Nanowire>.

## Setup

``` bash
 # Install
 git clone https://github.com/mikeyh30/nanowire-simulation.git
 cd nanowire-simulation
 pip install .

 # Setup configuration files
 cp config/globals-sample.yml globals.yml
 cp config/sim_parameters_units-sample.yml sim_parameters_units.yml
 cp config/sim_parameters_unitless-sample.yml sim_parameters_unitless.yml
```

Change the paths in ```globals.yml```  to point to the scratch, src, and hpc-control folder.
