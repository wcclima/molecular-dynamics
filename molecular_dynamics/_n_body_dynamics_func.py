
from typing import Callable, Tuple

import numpy as np

from ._compute_forces_func import pairwise_forces


def periodic_boundary_condition_xy(
        self, 
        x_coord: np.ndarray,
        y_coord: np.ndarray,
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Updates the molecules's Cartesian x- and 
    y-coordinates by imposing periodic boundary 
    conditions. 
    
    Returns the updated Cartesian coordinates array.
    
    Keyword arguments:
    
    x_coord (np.ndarray) -- the x Cartesian coordinate
    y_coord (np.ndarray) -- the y Cartesian coordinate
    """
    
    for coord in [x_coord, y_coord]:
        is_below = coord <= 0.
        is_above = coord >= self.box_width

        coord += is_below*self.box_width - is_above*self.box_width

    return x_coord, y_coord

def periodic_boundary_condition_x(
        self, 
        x_coord: np.ndarray,
        y_coord: np.ndarray,
        ) -> Tuple[np.ndarray, np.ndarray,]:
    """
    Updates the molecules's Cartesian x- and 
    y-coordinates by imposing periodic boundary 
    conditions. 
    
    Returns the updated Cartesian coordinates array.
    
    Keyword arguments:
    
    x_coord (np.ndarray) -- the x Cartesian coordinate
    y_coord (np.ndarray) -- the y Cartesian coordinate
    """
    
    is_below = x_coord <= 0.
    is_above = x_coord >= self.box_width

    x_coord += is_below*self.box_width - is_above*self.box_width

    return x_coord, y_coord

def diff_algorithm(
        coord: np.ndarray, 
        coord_prev: np.ndarray, 
        force: np.ndarray, 
        time_step: float, 
        mass: float = 1.
        ) -> np.ndarray:
    
    inv_mass = 1./mass
    time_step_squared = time_step*time_step
    
    return 2.*coord - coord_prev + inv_mass*force*time_step_squared

def verlet_method(
        self, 
        x_coord_prev: np.ndarray, 
        y_coord_prev: np.ndarray, 
        x_coord: np.ndarray, 
        y_coord: np.ndarray, 
        force_x: np.ndarray, 
        force_y: np.ndarray
        ) -> Tuple[float, float, float, float]:
    """
    Computes the positions and velocities of the particles 
    for the next time step using the Verlet, the past 
    positions and the total forces for each particle.
    
    Returns the positions and velocities of the particles 
    at the next time step.
    
    Keyword arguments:
    
    x_coord_prev (np.ndarray) -- x coordinate at the previous time step 
    y_coord_prev (np.ndarray) -- y coordinate at the previous time step
    x_coord (np.ndarray) -- x coordinate at the currrent time step 
    y_coord (np.ndarray) -- y coordinate at the currrent time step
    force_x (np.ndarray) -- force x component at the currrent time step
    force_y (np.ndarray) -- force y component at the currrent time step
    
    """

    x_coord_next = np.empty(self.n_molecules)
    y_coord_next = np.empty(self.n_molecules)
    
    x_coord_next[self.gas_molecules_index] = diff_algorithm(
        x_coord[self.gas_molecules_index], 
        x_coord_prev[self.gas_molecules_index], 
        force_x[self.gas_molecules_index], 
        self.time_step
        )
    y_coord_next[self.gas_molecules_index] = diff_algorithm(
        y_coord[self.gas_molecules_index], 
        y_coord_prev[self.gas_molecules_index], 
        force_y[self.gas_molecules_index], self.time_step
        )

    x_coord_next[self.brownian_particles_index] = diff_algorithm(
        x_coord[self.brownian_particles_index], 
        x_coord_prev[self.brownian_particles_index], 
        force_x[self.brownian_particles_index], 
        self.time_step,
        self.brownian_particle_mass
        )
    y_coord_next[self.brownian_particles_index] = diff_algorithm(
        y_coord[self.brownian_particles_index], 
        y_coord_prev[self.brownian_particles_index], 
        force_y[self.brownian_particles_index], 
        self.time_step,
        self.brownian_particle_mass        
        )

    vx = (x_coord_next - x_coord_prev)/(2.*self.time_step)
    vy = (y_coord_next - y_coord_prev)/(2.*self.time_step)
    
    return x_coord_next, y_coord_next, vx, vy
    

def time_evolution(self, 
                   x_coord_prev: np.ndarray, 
                   y_coord_prev: np.ndarray, 
                   x_coord: np.ndarray, 
                   y_coord: np.ndarray,
                   dist_func: Callable[[float, float], Tuple[float, float, float]],
                   wall_force_func: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the next step in the molecular dynamics 
    history assuming periodic boundary conditions in 
    the x- and y-directions.
    
    Returns the updated current positions and 
    velocities of the molecules in the gas.

    
    Keyword arguments:
    
    x_coord_prev (np.ndarray) -- x coordinate at the previous time step 
    y_coord_prev (np.ndarray) -- y coordinate at the previous time step
    x_coord (np.ndarray) -- x coordinate at the currrent time step 
    y_coord (np.ndarray) -- y coordinate at the currrent time step
    """
    
    force_x, force_y = (
        pairwise_forces(
            self,
            x_coord, 
            y_coord,
            dist_func
            )
        )
    
    force_y += - self.gravity_accel

    wall_force_x, wall_force_y = (
        wall_force_func(
            self, 
            x_coord, 
            y_coord
            )
        )

    force_x += wall_force_x
    force_y += wall_force_y 

    x_coord_next, y_coord_next, vx, vy = (
        verlet_method(
            self,
            x_coord_prev, 
            y_coord_prev, 
            x_coord, 
            y_coord, 
            force_x, 
            force_y
            )
        )
    
    x_coord_prev = x_coord
    x_coord = x_coord_next
    
    y_coord_prev = y_coord
    y_coord = y_coord_next
    
    return x_coord_prev, y_coord_prev, x_coord, y_coord, vx, vy


def dynamics(
        self,
        dist_func: Callable[[float, float], Tuple[float, float, float]],
        boundary_cond_func: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]],
        wall_force_func: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]
        ):


    new_data = np.empty((self.total_time_steps - 1, 4, self.n_molecules))
    self.data = np.concatenate((self.data, new_data))

    x_coord, y_coord, vx, vy = self.data[0]
        
    x_coord_prev, y_coord_prev = boundary_cond_func(
        self, 
        x_coord - vx*self.time_step,
        y_coord - vy*self.time_step
        )

    for time_step_i in range(1, self.total_time_steps):
        x_coord_prev, y_coord_prev, x_coord, y_coord, vx, vy = (
            time_evolution(
                self,
                x_coord_prev, 
                y_coord_prev, 
                x_coord, 
                y_coord,
                dist_func,
                wall_force_func
                )
            )
        
        x_coord, y_coord = boundary_cond_func(
            self, 
            x_coord, 
            y_coord
            )
        
        self.data[time_step_i] = np.array([x_coord, y_coord, vx, vy])