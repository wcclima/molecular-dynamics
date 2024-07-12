
from typing import Callable, Tuple, Union

import numpy as np

from ._n_body_dynamics_func import time_evolution


def thermostat(
        self, 
        new_temperature: Union[float, int], 
        extra_time_steps: int, 
        dist_func: Callable[[float, float], Tuple[float, float, float]], 
        boundary_cond_func: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]], 
        wall_force_func: Callable[[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray]]
        ):


    thermostat_factor = new_temperature/self.average_kT

    x_coord_prev, y_coord_prev = self.data[self.total_time_steps - 2][0:2]
    x_coord, y_coord, vx, vy= self.data[self.total_time_steps - 1]

    vx[self.gas_molecules_index] = thermostat_factor*vx[self.gas_molecules_index]
    vy[self.gas_molecules_index] = thermostat_factor*vy[self.gas_molecules_index]

    self.data[self.total_time_steps - 1][2] = vx
    self.data[self.total_time_steps - 1][3] = vy

    x_coord_prev[self.gas_molecules_index], y_coord_prev[self.gas_molecules_index] = (
        boundary_cond_func(
            self, 
            x_coord[self.gas_molecules_index] - (x_coord[self.gas_molecules_index] - x_coord_prev[self.gas_molecules_index])*thermostat_factor, 
            y_coord[self.gas_molecules_index] - (y_coord[self.gas_molecules_index] - y_coord_prev[self.gas_molecules_index])*thermostat_factor
            )
        )

    new_data = np.empty((extra_time_steps - 1, 4, self.n_molecules))

    for time_step_i in range(0,extra_time_steps - 1):
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
        
        x_coord, y_coord = (
            boundary_cond_func(
                self, 
                x_coord, 
                y_coord
                )
            )

        new_data[time_step_i] = np.array([x_coord, y_coord, vx, vy])

    self.total_time_steps += extra_time_steps - 1
    self.data = np.concatenate((self.data, new_data))