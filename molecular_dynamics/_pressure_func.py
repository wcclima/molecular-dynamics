
import numpy as np
import matplotlib.pyplot as plt

def measure_pressure(self):

    pressure = np.zeros(self.total_time_steps)

    for step_i in range(self.total_time_steps):
        x_coord = self.data[step_i][0][self.gas_molecules_index]
        vx = self.data[step_i][2][self.gas_molecules_index]

        dx_abs = np.abs(x_coord - 0.5*self.box_width)

        position_mask = dx_abs <= 0.25

        if position_mask.sum() != 0:
            pressure[step_i] = np.mean(np.abs(vx[position_mask]))/(self.time_step*self.box_width)

    self.mean_pressure = np.cumsum(pressure)/np.arange(1,self.total_time_steps + 1)

    print(f"Average pressure after {self.total_time_steps}: {self.mean_pressure.mean()}")

    fig = plt.figure(dpi = 150)
    plt.title("Pressure vs. time steps")
    plt.xlabel("time step")
    plt.ylabel("mean pressure")
    plt.plot(np.arange(self.total_time_steps), self.mean_pressure)
    plt.plot(np.arange(self.total_time_steps), pressure)
    plt.plot(np.arange(self.total_time_steps), np.array([self.mean_pressure.mean() for i in range(self.total_time_steps)]), 'k:')

    