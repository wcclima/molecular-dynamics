
from typeguard import typechecked
from typing import Union

import numpy as np

from ._initial_cond_func import generate_initial_conditions
from ._plot_func import plot_initial_conditions_bhw, plot_initial_forces_bhw
from ._animation_func import gas_animation_bhw , vel_distribution_animation
from ._stats_func import vel_distribution, measure_kT
from ._n_body_dynamics_func import dynamics
from ._thermostat_func import thermostat
from ._compute_forces_func import shortest_distance_x_pbc, bottom_wall_force
from ._n_body_dynamics_func import periodic_boundary_condition_x


class MissingPrecedingMethodCallError(Exception):
    def __init__(self, message):
        super().__init__(message)


__all__ = ["MolecularDynamicsBHW"]


class MolecularDynamicsBHW(object):

    @typechecked
    def __init__(
        self, 
        box_width: Union[float, int], 
        gas_density: Union[float, int], 
        time_step: Union[float, int], 
        initial_temperature: Union[float, int],
        gravity_accel: Union[float, int]
        ):


        self.box_width = box_width
        self.gas_density = gas_density
        self.time_step = time_step
        self.initial_temperature = initial_temperature
        self.gravity_accel = gravity_accel

        self.cell_size = np.sqrt(1./self.gas_density)
        self.n_cells = int(np.sqrt(self.gas_density)*self.box_width)
        self.n_molecules = self.n_cells*self.n_cells
        self.n_gas_molecules = self.n_molecules
        self.n_brownian_particles = 0
        self.brownian_particle_mass = 1.
        self.brownian_particle_size = None
        
        self.initial_state_flag = False
        self.dynamics_flag = False
        self.vel_dist_flag = False
        self.measured_kT_flag = False
        self.track_brownian_particles_flag = False

        self.gas_state = None
        self.brownian_particles_state = None
        self.brownian_motion_track = None
        self.vel_hist = None

    
    def display_parameters(self):
        """
        Prints all the relevant parameters in the model.
        """
        
        print("Model's parameters:")
        print("---"*12)
        print(f"Box width: {self.box_width}")
        print(f"Gas density: {self.gas_density}")
        print(f"Time step: {self.time_step}")
        print(f"Initial temperature: {self.initial_temperature}")
        print(f"Gravity acceleration: {self.gravity_accel}")
        print(f"Average molecule separation: {np.round(self.cell_size, 2)}")
        print(f"Number of molecules: {self.n_molecules}")


    @typechecked
    def generate_initial_conditions(
        self, 
        random_initial_position: bool = True
        ):

        self.initial_state_flag = True

        return generate_initial_conditions(
            self, 
            random_initial_position
            )
    

    def plot_initial_conditions(self):

        if not self.initial_state_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_initial_conditions' must be called before 'plot_initial_conditions' is called."
            )
        
        return plot_initial_conditions_bhw(self)


    def plot_initial_forces(self):


        if not self.initial_state_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_initial_conditions' must be called before 'plot_initial_forces' is called."
            )
        
        return plot_initial_forces_bhw(self)
    

    @typechecked
    def generate_dynamics(
        self, 
        total_time_steps: int
        ):

        if not self.initial_state_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_initial_conditions' must be called before 'generate_dynamics' is called."
            )

        self.total_time_steps = total_time_steps
        self.dynamics_flag = True

        return dynamics(
            self, 
            shortest_distance_x_pbc,
            periodic_boundary_condition_x,
            bottom_wall_force
            )


    @typechecked
    def generate_gas_animation(
        self, 
        exit_file_name: str = "molecular_dynamics.mp4", 
        fps: Union[float, int] = 50, 
        dpi: Union[float, int] = 150
        ):

        if not self.dynamics_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_dynamics' must be called before 'generate_gas_animation' is called."
            )
        
        return gas_animation_bhw(
            self, 
            exit_file_name = exit_file_name, 
            fps = fps, 
            dpi = dpi
            )


    @typechecked
    def generate_velocity_distribution(
            self, 
            bins: int = 500, 
            bin_size: Union[float, int] = 1e-2
            ):

        if not self.dynamics_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_dynamics' must be called before 'generate_velocity_distribution' is called."
            )

        self.vel_dist_flag = True

        return vel_distribution(
            self, 
            bins, 
            bin_size
            )


    @typechecked
    def generate_velocity_distribution_animation(
            self, 
            exit_file_name: str = "velocity_distribution.mp4"
            ):
        
        if not self.vel_dist_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_velocity_distribution' must be called before 'generate_velocity_distribution_animation' is called."
            )

        return vel_distribution_animation(
            self, 
            exit_file_name
            )        
    

    @typechecked
    def measure_kT(
        self, 
        method: str = "instantaneous fitting"
        ):

        if not self.vel_dist_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_velocity_distribution' must be called before 'measure_kT' is called."
            )
        
        self.measured_kT_flag = True
        
        return measure_kT(
            self, 
            method
            )
    
    
    def change_temperature(
            self, 
            new_temperature: Union[float, int], 
            extra_time_steps: int
            ):

        if not self.dynamics_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_dynamics' must be called before 'change_temperature' is called."
            )
        
        if not self.measured_kT_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'measure_kT' must be called before 'change_temperature' is called."
            )
        
        return thermostat(
            self, 
            new_temperature, 
            extra_time_steps,
            shortest_distance_x_pbc, 
            periodic_boundary_condition_x, 
            bottom_wall_force
            )
