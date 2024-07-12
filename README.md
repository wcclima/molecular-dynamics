# Molecular Dynamics
Repo for the molecular dynamics project.

## 1 - Objective

The goal of this project is to simulate a gas by modelling its molecules as classical particles interacting with each other via the Lennard-Jones potential. The simulation allows us to produce animations depicting the gas molecules's movement, to see the gas achieve its thermal equilibrium, to measure its temperature and to change its temperature to show phase transition. 

## 2 - Repo organisation

**`notebooks/`: Notebooks demonstrating the simulations**
- `molecular_dynamics_notebook.ipynb` : Notebook explaining the basic details of the simulation
- `brownian_motion_notebook.ipynb` : Notebook discussing the Brownian motion implementation
- `bottom_hard_wall_notebook.ipynb` : Notebook discussing the implementation of the bottom-hard-wall boundary condition

**`molecular_dynamics/`: Molecular dynamics module**

See module architecture.


## 3 - Module architecture

Description of the module `molecular_dynamics` architecture.

- `molecular_dynamics/__init__.py`
  - Initialises the module.
  - Imports the MolecularDynamicsPBC class, with periodic boundary condition (PBC).
  - Imports the MolecularDynamicsBHW class, with bottom hard wall boundary condition (BHW).

- `molecular_dynamics/molecular_dynamics_pbc.py`: defines the MolecularDynamicsPBC class with the methods
  -  `display_parameters`
  -  `generate_initial_conditions`, requires functions from `_initial_cond.py`; 
  -  `plot_initial_conditions`, requires functions from `_plot_methods.py`;
  -  `plot_initial_forces`, requires functions from `_plot_methods.py`;
  -  `generates_dynamics`, requires functions from `_n_body_dynamics.py`;
  -  `generate_gas_animation`, requires functions from `_animation_methods.py`;
  -  `generate_velocity_distribution`, requires functions from `_stats_methods.py`;
  -  `generate_velocity_distribution_animation`, requires functions from `_animation_methods.py`;
  -  `measure_kT`, requires functions from `_stats_methods.py`;
  -  `change_temperature`, requires functions from `_thermostat.py`.
- `molecular_dynamics/molecular_dynamics_bhw.py`: defines the MolecularDynamicsBHW class with the same methods as above.
- `molecular_dynamics/_animation_func.py`: functions producing the animations of the gas moviment and the velocity histogram, used in the `generate_gas_animation` and `generate_velocity_distribution_animation` methods.
- `molecular_dynamics/_compute_forces_func.py`: functions computing the molecules pairwise total forces and the bottom hard wall forces, used in the `plot_initial_forces` method and in `_n_body_dynamics.py` functions.
- `molecular_dynamics/_initial_cond_func.py`: function generating the molecules's initial conditions, used in the `generate_initial_conditions` method. 
- `molecular_dynamics/_n_body_dynamics_func.py`: functions generating the gas dynamics, used in the `generates_dynamics` and `change_temperature` methods.
- `molecular_dynamics/_plot_func.py`: functions generating the plots, used in the `plot_initial_conditions` and `plot_initial_forces` methods.
- `molecular_dynamics/_stats_func.py`: functions implementing the temperature measurement by generating the velocity distribution and fitting it with the Maxwell-Boltzmann distribution or using the equipartition theorem, used in the `measure_kT` and `generate_velocity_distribution` methods.
- `molecular_dynamics/_thermostat_func.py`: functions that instantaneously change the molecules's velocities, thus changing the gas temperature after interation, used in the `change_temperature` method.

## 4 - Results

## 5 - References

1. N. J. Giordano and H Nakanishi, *Computational Physics* (Pearson Prentice Hall, New Jersey, 2006).
2. R.H. Landau, M.J. Páez and C.C. Bordeianu, *Computation Physics: Problem Solving with Python* (Wiley-VCH, Weinheim, 2015). 
