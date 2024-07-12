
from typing import Callable, Tuple

import numpy as np

def force_func(a: float, 
               b: float, 
               epsilon: float, 
               x: float
               ) -> float:
    """
    Returns the modulus of the force 
    divide by the Euclidean distance.

    Keyword arguments:

    a (float) -- parameter of the repulsive term in the Lennar-Jones potential
    b (float) -- parameter of the attractive term in the Lennar-Jones potential
    epsilon (float) -- coupling parameter
    x (float) -- inverse distance squared    
    """
    return 24.*epsilon*(2.*a*x**7 - b*x**4)

def shortest_distance_xy_pbc(
        self, 
        x_coord_distance: float, 
        y_coord_distance: float
        ) -> Tuple[float, float, float]:
    """
    Computes the shortest coordinate distance between
    two particles using periodic boundary conditions 
    in the x- and y-directions.

    Returns the updated Cartesian coordinate distances in 
    the x- and y-directions and the Euclidean squared 
    distance.
    
    Keyword arguments:
    
    x_coord_distance (float) -- Cartesian x-coordinate distance.
    y_coord_distance (float) -- Cartesian y-coordinate distance. 
    """
    
    if (np.abs(x_coord_distance) > 0.5*self.box_width):
        x_coord_distance = x_coord_distance - np.copysign(self.box_width, x_coord_distance)

    if (np.abs(y_coord_distance) > 0.5*self.box_width):
        y_coord_distance = y_coord_distance - np.copysign(self.box_width, y_coord_distance)
        
    squared_distance = x_coord_distance**2 + y_coord_distance**2

    return x_coord_distance, y_coord_distance, squared_distance


def shortest_distance_x_pbc(
        self, 
        x_coord_distance: float, 
        y_coord_distance: float
        ) -> Tuple[float, float, float]:
    """
    Computes the shortest coordinate distance between
    two particles using periodic boundary conditions 
    in the x-direction.

    Returns the updated Cartesian coordinate distances in 
    the x- and y-directions and the Euclidean squared 
    distance.
    
    Keyword arguments:
    
    x_coord_distance (float) -- Cartesian x-coordinate distance.
    y_coord_distance (float) -- Cartesian y-coordinate distance. 
    """
    
    if (np.abs(x_coord_distance) > 0.5*self.box_width):
        x_coord_distance = x_coord_distance - np.copysign(self.box_width, x_coord_distance)
        
    squared_distance = x_coord_distance**2 + y_coord_distance**2

    return x_coord_distance, y_coord_distance, squared_distance


def zero_force(
        self, 
        x_coord: np.ndarray, 
        y_coord: np.ndarray
        ) -> Tuple[np.ndarray, np.ndarray]:

    return np.zeros(self.n_molecules), np.zeros(self.n_molecules)

def bottom_wall_force(
        self, 
        x_coord: np.ndarray, 
        y_coord: np.ndarray
        ) -> Tuple[np.ndarray, np.ndarray]:
    
    x_wall_cells = np.array([i*self.box_width/50 for i in range(50)])
    force_wall_x = np.zeros(self.n_molecules)
    force_wall_y = np.zeros(self.n_molecules)
    
    for i in range(self.n_molecules):
        for j in range(50):

            dx, dy, squared_dist = shortest_distance_x_pbc(
                self, 
                x_coord[i] - x_wall_cells[j],
                y_coord[i]
                )            
            
            if (squared_dist < 0.25):

                inv_squared_dist = 1./squared_dist
                force_wall_x[i] += 1200.*(0.09*inv_squared_dist)**7*dx
                force_wall_y[i] += 1200.*(0.09*inv_squared_dist)**7*dy
                
                
    return force_wall_x, force_wall_y


def pairwise_forces(
        self, 
        x_coord: np.ndarray, 
        y_coord: np.ndarray, 
        dist_func: Callable[[float, float], Tuple[float, float, float]]
        ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the forces and potential energy between any pair 
    of particles using the Lennard-Jonnes potential and the 
    shortest distance as defined by the argument 'dist_func'.
    
    Returns the x and y components of the force between the 
    particles.
    
    Keyword arguments:
    
    x_coord (np.ndarray) -- x coordinate at the currrent time step 
    y_coord (np.ndarray) -- y coordinate at the currrent time step
    dist_func
    """
    
    
    force_x = np.zeros(self.n_molecules)
    force_y = np.zeros(self.n_molecules)

    for i in range(self.n_gas_molecules - 1):
        gas_molecule_i = self.gas_molecules_index[i]
        for j in range(i + 1, self.n_gas_molecules):
            gas_molecule_j = self.gas_molecules_index[j]

            dx, dy, dist_squared = dist_func(
                self, 
                x_coord[gas_molecule_i] - x_coord[gas_molecule_j],
                y_coord[gas_molecule_i] - y_coord[gas_molecule_j]
                )
                
            if (dist_squared < 9.): #avoids computing negligible forces
                if (dist_squared < 1e-3): #avoids 0 denominator
                    dist_squared = 1e-3

                inv_dist_squared = 1./dist_squared
                f = force_func(1., 1., 1., inv_dist_squared)
            else:
                f = 0.
                                  
            force_x[gas_molecule_i] = f*dx + force_x[gas_molecule_i]
            force_x[gas_molecule_j] = -f*dx + force_x[gas_molecule_j]
                
            force_y[gas_molecule_i] = f*dy + force_y[gas_molecule_i] 
            force_y[gas_molecule_j] = -f*dy + force_y[gas_molecule_j]

    for i in range(self.n_brownian_particles - 1):
        brownian_particle_i = self.brownian_particles_index[i]
        for j in range(i + 1, self.n_brownian_particles):
            brownian_particle_j = self.brownian_particles_index[j]

            dx, dy, dist_squared = dist_func(
                self, 
                x_coord[brownian_particle_i] - x_coord[brownian_particle_j],
                y_coord[brownian_particle_i] - y_coord[brownian_particle_j]
                )
                
            if (dist_squared < 9.): #avoids computing negligible forces
                if (dist_squared < 1e-3): #avoids 0 denominator
                    dist_squared = 1e-3

                inv_dist_squared = 1./dist_squared
                f = force_func((2.*self.brownian_particle_size)**6,1.,1e1,inv_dist_squared)
            else:
                f = 0.
                                  
            force_x[brownian_particle_i] = f*dx + force_x[brownian_particle_i]
            force_x[brownian_particle_j] = -f*dx + force_x[brownian_particle_j]
                
            force_y[brownian_particle_i] = f*dy + force_y[brownian_particle_i] 
            force_y[brownian_particle_j] = -f*dy + force_y[brownian_particle_j]

    for i in range(self.n_brownian_particles):
        brownian_particle_i = self.brownian_particles_index[i]
        for j in range(self.n_gas_molecules):
            gas_molecule_j = self.gas_molecules_index[j]

            dx, dy, dist_squared = dist_func(
                self, 
                x_coord[brownian_particle_i] - x_coord[gas_molecule_j],
                y_coord[brownian_particle_i] - y_coord[gas_molecule_j]
                )
                
            if (dist_squared < 9.): #avoids computing negligible forces
                if (dist_squared < 1e-3): #avoids 0 denominator
                    dist_squared = 1e-3

                inv_dist_squared = 1./dist_squared
                f = force_func(self.brownian_particle_size**6*1.,1.,1e1,inv_dist_squared)
            else:
                f = 0.
                                  
            force_x[brownian_particle_i] = f*dx + force_x[brownian_particle_i]
            force_x[gas_molecule_j] = -f*dx + force_x[gas_molecule_j]
                
            force_y[brownian_particle_i] = f*dy + force_y[brownian_particle_i] 
            force_y[gas_molecule_j] = -f*dy + force_y[gas_molecule_j]
            
            
    return force_x, force_y