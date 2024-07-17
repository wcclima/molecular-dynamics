
import numpy as np

import matplotlib.pyplot as plt

from ._compute_forces_func import pairwise_forces, shortest_distance_xy_pbc, shortest_distance_x_pbc

def plot_initial_conditions_pbc(self):

    x_0, y_0, vx_0, vy_0 = self.data[0]

    fig = plt.figure(dpi=150)

    plt.axis([0., self.box_width, 0., self.box_width])

    plt.plot(x_0[self.brownian_particles_index], y_0[self.brownian_particles_index], 
            marker = "o", 
            linestyle = " ",
            markerfacecolor = "orange", 
            markeredgecolor = "orange",  
            markersize=10,
            zorder = 1
        )

    plt.plot(x_0[self.gas_molecules_index], y_0[self.gas_molecules_index], 
            marker = "o", 
            linestyle = " ", 
            markerfacecolor = "cornflowerblue", 
            markeredgecolor = "cornflowerblue", 
            markersize=5,
            zorder = 1
            )
    for i in range(self.n_molecules):
        plt.arrow(
            x_0[i],
            y_0[i],
            vx_0[i],
            vy_0[i],
            length_includes_head = True,
            head_width = 0.1, 
            color = "black",
            zorder = 2
            )

    plt.title("Initial positions and velocities of the molecules in the gas")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.show()


def plot_initial_conditions_bhw(self):

    x_0, y_0, vx_0, vy_0 = self.data[0]

    x_wall = np.array([i*self.box_width/10. for i in range(11)])

    fig = plt.figure(dpi=150)

    plt.axis([0., self.box_width, -2., 1.5*self.box_width])

    plt.plot(x_wall, 0.*x_wall, "k")
    for i in range(-1,41):
        x_aux = np.array([j*self.box_width/40. + i*self.box_width/40. for j in range(3)])
        y_aux = 2.*x_aux - 2.*(0.5*i + 1)
        plt.plot(x_aux, y_aux, ":k")

    plt.plot(x_0, y_0, 
            marker = "o", 
            linestyle = " ", 
            markerfacecolor = "cornflowerblue", 
            markeredgecolor = "cornflowerblue", 
            markersize=5,
            zorder = 1
            )
    for i in range(self.n_molecules):
        plt.arrow(
            x_0[i],
            y_0[i],
            vx_0[i],
            vy_0[i],
            length_includes_head = True,
            head_width = 0.1, 
            color = "black",
            zorder = 2
            )

    plt.title("Initial positions and velocities of the molecules in the gas")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.show()


def plot_initial_forces_pbc(self):

    x_0, y_0 = self.data[0][0:2]

    fx, fy = pairwise_forces(
        self, 
        x_0, 
        y_0,
        shortest_distance_xy_pbc
        )

    fig = plt.figure(dpi=150)

    plt.axis([0., self.box_width, 0., self.box_width])

    plt.plot(x_0[self.brownian_particles_index], y_0[self.brownian_particles_index], 
            marker = "o", 
            linestyle = " ",
            markerfacecolor = "orange", 
            markeredgecolor = "orange",  
            markersize=10,
            zorder = 1
        )

    plt.plot(x_0[self.gas_molecules_index], y_0[self.gas_molecules_index], 
            marker = "o", 
            linestyle = " ", 
            markerfacecolor = "cornflowerblue", 
            markeredgecolor = "cornflowerblue", 
            markersize=5,
            zorder = 1
            )
    for i in range(self.n_molecules):
        plt.arrow(
            x_0[i],
            y_0[i],
            fx[i],
            fy[i],
            length_includes_head = True, 
            head_width = 0.1, 
            color='black'
            )

    plt.title("Initial position and total force of the molecules in the gas")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.show()


def plot_initial_forces_bhw(self):

    x_0, y_0 = self.data[0][0:2]
    x_wall = np.array([i*self.box_width/10. for i in range(11)])

    fx, fy = pairwise_forces(
        self, 
        x_0, 
        y_0,
        shortest_distance_x_pbc
        )
    
    fy += - self.gravity_accel

    fig = plt.figure(dpi=150)

    plt.axis([0., self.box_width, -2., 1.5*self.box_width])

    plt.plot(x_wall, 0.*x_wall, "k")
    for i in range(-1,41):
        x_aux = np.array([j*self.box_width/40. + i*self.box_width/40. for j in range(3)])
        y_aux = 2.*x_aux - 2.*(0.5*i + 1)
        plt.plot(x_aux, y_aux, ":k")

    plt.plot(x_0, y_0, 
            marker = "o", 
            linestyle = " ", 
            markerfacecolor = "cornflowerblue", 
            markeredgecolor = "cornflowerblue", 
            markersize=5,
            zorder = 1
            )
    for i in range(self.n_molecules):
        plt.arrow(
            x_0[i],
            y_0[i],
            fx[i],
            fy[i],
            length_includes_head = True, 
            head_width = 0.1, 
            color='black'
            )

    plt.title("Initial position and total force of the molecules in the gas")
    plt.xlabel("$x$")
    plt.ylabel("$y$")

    plt.show()
