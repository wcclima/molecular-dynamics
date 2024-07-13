
from typing import Union

import numpy as np
from numpy import random

def gaussian_random_array(
        shape: Union[int, tuple],
        trials: int,
        ) -> np.ndarray:

    """
    Returns a numpy array that is an approximate random 
    Gaussian variable with zero mean and elements in the 
    range [-0.5, 0.5].

    Keyword arguments:

    shape (int, tuple) -- the shape of the array
    trials (int) -- length of the random sequence
    
    """
    s = np.zeros(shape)
    
    for i in range(trials):
        s += np.random.rand(shape)
    
    return s/trials - 0.5

def generate_initial_conditions(
        self, 
        random_initial_position: bool = True
        ) -> None:
    
    """
    Generates a random initial condition for the position and velocity 
    of the molecules. The positions can be either on the vertice of a
    square lattice or randomly displaced from that position in the xy 
    directions by a*d by a quarter of the lattice cell size. 
    
    Returns the initial xy molecules positions and their xy velocity 
    components.

    Keyword arguments:

    random_initial_position (bool) -- True is the standard value. It 
    determines whether the lattice is regular (False) or randomized (True).

    """
    
    self.gas_molecules_index = [i for i in range(self.n_molecules)]
    self.brownian_particles_index = []
    for i in range(self.n_brownian_particles):
        bp_index = random.randint(self.n_molecules)
        self.gas_molecules_index.remove(bp_index)
        self.brownian_particles_index.append(bp_index)

    vx = 12.*gaussian_random_array(self.n_molecules, 12)*np.sqrt(self.initial_temperature)
    vy = 12.*gaussian_random_array(self.n_molecules, 12)*np.sqrt(self.initial_temperature)
    
    cell_array = np.linspace(
        0.5*self.cell_size, 
        self.box_width - 0.5*self.cell_size, 
        self.n_cells
        )

    x_coord, y_coord = (
        np.meshgrid(cell_array, cell_array)
        + random_initial_position*(
            np.random.rand(2, 
                            self.n_cells, 
                            self.n_cells
                            ) - 0.5
                            )*0.5*self.cell_size
        )
    
    x_coord = x_coord.flatten()
    y_coord = y_coord.flatten()

    self.data = np.empty((1, 4, self.n_molecules))
    self.data[0] = np.array([x_coord, y_coord, vx, vy])
