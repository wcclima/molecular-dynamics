
from typeguard import typechecked
from typing import Optional

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

    """
    Simulates the molecular dynamics of a gas with bottom hard
    wall (BHW) conditions.

    MolecularDynamicsBHW simulates the dynamic of the molecules 
    of a bi-dimensional gas interacting via a 6-12 Lennard-Jones 
    potential and assuming periodic boundary conditions in the 
    x-direction and hard walls for y = 0 and open boundary condition
    for y > 0. The hard wall is modelled by an array of static 
    particles that interact with the gas molecules with the repulsive 
    12 part of the Lennard-Jones potential. The model also allows for
    a constant homogenous force field in the negative y-direction
    to simulate the gravitational acceleration. The dynamics is 
    generated with Newton's 2nd Law implemented using Verlet's 
    method. All the parameters scaled by the gas molecule's
    mass, size and coupliong strength. 

    Attributes:
        n_molecules (int): 
            The total number of molecules.

        data (np.ndarray):
            The array containing the state history of all
            molecules in the system over time.

        total_time_steps (int):
            Total number of time steps recorded.

        average_KT (float):
            The current equilibrium temperature in energy units.

    Methods:
        display_parameters()
            Prints all the relevant parameters in the model.

        generate_initial_conditions(random_initial_position):
            Generates the initial positions and velocities of the 
            molecules in the system.

        plot_initial_conditions():
            Plots the initial positions and velocities of the 
            molecules in the system.

        plot_initial_forces():
            Plots the initial positions and forces of the 
            molecules in the system.

        generate_dynamics(total_time_steps):
            Generates the dynamics of the molecules in the system.
            The dynamics is generated according to Newton's 2dn Law
            implemented via Verlet's method and assuming a 6-12
            Lennard-Jones potential, periodic boundary condition
            in the x-direction and botton hard wall condition in 
            the y-direction.

        generate_gas_animation(exit_file_name, fps, dpi):
            Generates the animation of the molecule's positions. 
            It saves the animation in the mp4 format.

        generate_velocity_distribution(bins, bin_size):
            Generates the histogram corresponding to the velocity 
            modulus distribution of the molecules for each time 
            step. 

        generate_velocity_distribution_animation(exit_file_name, average_window):
            Generates the animation of the velocity modulus histograms
            for the molecules. It saves the animation in the mp4 format.

        measure_kT(method):
            Measures the gas temperature by either fitting the velocity 
            distribution with the Maxwell-Boltzmann distribution for the
            velocity modulus or by computing the kinectic energy average
            and usinging the equipartition theorem.

        measure_pressure():
            Measures the pressure exerted by the molecules on an
            imaginary perfectly reflective wall.

        change_temperature():
            Changes the molecules temperature by rescalling the 
            molecule's velocity and evolving the system until the new 
            equilibrium is reached.

    """

    @typechecked
    def __init__(
        self, 
        box_width: float = 20., 
        gas_density: float = 0.25, 
        time_step: float = 5e-3, 
        initial_temperature: float = 1.,
        gravity_accel: float = 0.5
        ):
        """
        Initialises MolecularDynamicsBHW with parameters.

        Keyword arguments:
            box_width (float, default = 20.):
                The length of the bi-dimensional volume side.

            gas_density (float, default = 0.25):
                The gas density.

            time_step (float, default = 5e-3):
                The size of the time step.

            initial_temperature (float, default = 1.):
                The gas initial temperature kT, where k is Boltzmann 
                constant. It is the variance of the x- and y-components 
                of the initial velocities.

            gravity_accel (float, default = 0.5):
                The gravitational acceleration. 
        """

        # class parameters
        self.box_width = box_width
        self.gas_density = gas_density
        self.time_step = time_step
        self.initial_temperature = initial_temperature
        self.gravity_accel = gravity_accel

        #parameters that will be inputed after the class is initialised
        self.cell_size = np.sqrt(1./self.gas_density)
        self.n_cells = int(np.sqrt(self.gas_density)*self.box_width)
        self.n_molecules = self.n_cells*self.n_cells
        self.n_gas_molecules = self.n_molecules
        self.n_brownian_particles = 0
        self.brownian_particle_mass = 1.
        self.brownian_particle_size = None
        
        #flags to check if the methods have been called in the right order
        self._initial_state_flag = False
        self._dynamics_flag = False
        self._vel_dist_flag = False
        self._measured_kT_flag = False
        self._track_brownian_particles_flag = None

        #variables to keep the animation frames
        self._vel_bins = None
        self._vel_frac_count = None
        self._gas_state = None
        self._brownian_particles_state = None
        self._brownian_motion_track = None
        self._vel_hist = None

    
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
        """
        Generates the initial conditions of the molecules. 
        The positions are set as the vertices of a square 
        grid with cell size side as the avera particle separation. 
        The x- and y-components of the velocities are generated 
        as Gaussian random variables with zero mean and 
        variance given by the initial temperature.

        A random displacement is added to the positions 
        on square lattice if random_initial_position is 
        True.

        Keyword arguments:
            random_initial_position (bool, default = True):
                If True, the initial positions is made random
                by adding a random displacement to the 
                positions at the square grid.            

        """

        self.initial_state_flag = True

        return generate_initial_conditions(
            self, 
            random_initial_position
            )
    

    def plot_initial_conditions(
            self
            ):
        """
        Plots the initial positions and velocities of the 
        molecules in system. The velocities are plotted as 
        black arrows sitting on the molecules.  
        """

        if not self.initial_state_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_initial_conditions' must be called before 'plot_initial_conditions' is called."
            )
        
        return plot_initial_conditions_bhw(self)


    def plot_initial_forces(
            self
            ):
        """
        Plots the initial positions and forces of the 
        molecules. The forces are plotted as black arrows 
        sitting on the molecules.  
        """

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
        """
        Generates the gas molecules dynamics according to 
        Newton's 2nd Law implemented using Verlet's 
        method and the Lennard-Jones potential.

        Keyword arguments:
            total_time_steps (int):
                Total number of total steps.
        """

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
        fps: int = 50, 
        dpi: int = 150
        ):
        """
        Generates the animation of the molecule's positions. 
        It saves the animation in the mp4 format.

        Keyword arguments:
            exit_file_name (str, default = "molecular_dynamics.mp4"):
                Name of the exit file.

            fps (int, default = 50):
                Frames per second.

            dpi (int, default = 150):
                Dots per inch.            

        """

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
            bin_size: float = 1e-2
            ):
        """
        Generates the velocity modulus distribution as 
        a histogram for each time step.

        Keyword arguments:
            bins (int, default = 500):
                Number of bins in the histogram.

            bin_size (float, default = 1e-2):
                Size of the individual bins in the histogram.
        """


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
            exit_file_name: str = "velocity_distribution.mp4",
            average_window: int = None
            ):
        
        if average_window is None:
            average_window = self.total_time_steps

        """
        Generates the animation of the velocity modulus 
        distribution of the gas molecules as a histogram.
        The frames generated from the velocity modulus 
        distribution histogram histograms averaged over time,
        up to a certain number of past time steps set by 
        average_widow. It saves the animation in the mp4 format.

        Keyword arguments:
            exit_file_name (str, default = "velocity_distribution.mp4"):
                Name of the exit file.
            average_window (int, default = total_time_steps):
                The temporal moving average window.
        """
        
        if not self.vel_dist_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_velocity_distribution' must be called before 'generate_velocity_distribution_animation' is called."
            )

        return vel_distribution_animation(
            self, 
            exit_file_name,
            average_window
            )        
    

    @typechecked
    def measure_kT(
        self, 
        method: str = "equipartition",
        average_window: Optional[int] = None
        ):
        """
        Measures the gas molecules temperature in energy 
        scale. It prints the value for kT at the last
        recorded time step.
        
        The temperature can be measured by three different
        methods.
        - instantaneous fitting: Fits the velocity modulus 
        distribution at the last time step with the 
        Maxwell-Boltzmann distribution after a log-log
        transformation. The method returns the kT as the 
        average values between the constant and slope 
        regression coefficients and the plots of the velocity 
        distribution and the corresponding fitted curve. 
        
        - averaged fitting: Fits the time averaged velocity 
        modulus distribution at the last time step with the 
        Maxwell-Boltzmann distribution after a log-log
        transformation. The time average is perfomed by averaging 
        the velocity histograms for up to a certain time window 
        set by average_window parameter and ending at the 
        current time step. The method returns the kT as the 
        average values between the constant and slope 
        regression coefficients the plots of the velocity 
        distribution and the corresponding fitted curve.

        - equipartition: Computes the average kinetic energy 
        at the last time step and used the equipartition 
        theorem in two dimensions to obtain kT.

        Keyword arguments:
            method (str, default = 'equipartition'):
                The method used to measure kT.

            average_window (int, optional, default = None):
                The temporal average window size. Only used with
                the 'averaged fitting' method.
        """

        if not self.vel_dist_flag:
            raise MissingPrecedingMethodCallError(
                "Error: 'generate_velocity_distribution' must be called before 'measure_kT' is called."
            )
        
        self.measured_kT_flag = True
        
        return measure_kT(
            self, 
            method,
            average_window
            )
    
    
    def change_temperature(
            self, 
            new_temperature: float, 
            extra_time_steps: int
            ):
        """
        Changes the temperature of the gas molecules 
        by rescaling their velocities by a factor 
        sqrt(new_temperature/old_temperature) and then 
        evolving the system for an extra number of time steps.

        Keyword arguments:
            new_temperature (float):
                The new temperature the gas should approach 
                in energy units.

             extra_time_steps (int):
                The additional time steps to evolve the system
                until it reaches the new equilibrium. 
        """

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
