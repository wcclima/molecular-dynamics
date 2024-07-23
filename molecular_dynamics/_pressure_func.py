
import numpy as np
import matplotlib.pyplot as plt

from ._stats_func import time_average_func


def measure_pressure(
        self, 
        average_window: int
        ):
    
    """
    Computes the pressure exerted by the molecules on a virtual perfectly 
    refletive wall.

    
    Keyword arguments:

        average_window (int):
            The time average window in time steps.
    """

    pressure = np.zeros(self.total_time_steps)

    for step_i in range(self.total_time_steps):
        x_coord = self.data[step_i][0][self.gas_molecules_index]
        vx = self.data[step_i][2][self.gas_molecules_index]

        dx_abs = np.abs(x_coord - 0.5*self.box_width)

        position_mask = dx_abs <= 0.25

        if position_mask.sum() != 0:
            pressure[step_i] = np.mean(np.abs(vx[position_mask]))/(self.time_step*self.box_width)


    self.pressure_averaged = np.empty(self.total_time_steps)

    for step_i in range(self.total_time_steps):
        self.pressure_averaged[step_i] = time_average_func(pressure, average_window, step_i)

    print(f"Average pressure after {self.total_time_steps}: {self.pressure_averaged[-1]}")

    fig = plt.figure(dpi = 150)
    plt.title("Pressure vs. time steps")
    plt.xlabel("time step")
    plt.ylabel("mean pressure")
    plt.plot(np.arange(self.total_time_steps), self.pressure_averaged)
    
