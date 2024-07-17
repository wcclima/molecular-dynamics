
import numpy as np

from typing import Tuple

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure

import functools

import warnings


def init_state(self) -> Figure:
    """
    Set the gas's initial state.
    
    Returns a figure of the scatter plot of 
    the molecules's initial positions.
    """

    x_coord, y_coord = self.data[0][0:2]
    self._gas_state.set_data([x_coord[self.gas_molecules_index], y_coord[self.gas_molecules_index]])
    self._brownian_particles_state.set_data([x_coord[self.brownian_particles_index], y_coord[self.brownian_particles_index]])
    
    return self._gas_state, self._brownian_particles_state


def mask_periodic_plot(
        x_coord: np.ndarray, 
        y_coord: np.ndarray
        ) -> Tuple[np.ndarray, np.ndarray]:
    
    """
    Breaks a trajectory in periodic boundary conditions 
    by masking its jumps.

    Returns the masked x and y coordinates arrays.

    
    Keyword arguments:
        
    x_coord (np.ndarray) -- the trajectory's x-coordinates
    y_coord (np.ndarray) -- the trajectory's y-coordinates
    """
    
    abs_dx = np.abs(np.diff(x_coord))
    abs_dy = np.abs(np.diff(y_coord))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mask_x_coord = np.hstack([abs_dx > np.nanmean(abs_dx) + 3.*np.nanstd(abs_dx), [False]])
        mask_y_coord = np.hstack([abs_dy > np.nanmean(abs_dy) + 3.*np.nanstd(abs_dy), [False]])
    
    x_coord = np.ma.MaskedArray(x_coord, mask_x_coord)
    y_coord = np.ma.MaskedArray(y_coord, mask_y_coord)
    
    return x_coord, y_coord

def animate_gas(
        self, 
        frame_i: int
        ) -> Figure:

    """
    Set the gas's state at the i-th time step.
    
    Returns a figure of the scatter plot of 
    the molecules's positions at the i-th step.


    Keyword arguments:
        
    frame_i (int) -- the index of the i-th frame 
    """

    x_coord, y_coord = self.data[frame_i][0:2]
    self._gas_state.set_data(x_coord[self.gas_molecules_index], y_coord[self.gas_molecules_index])
    self._brownian_particles_state.set_data([x_coord[self.brownian_particles_index], y_coord[self.brownian_particles_index]])

    if self._track_brownian_particles_flag:
        for brownian_particle_i in range(self.n_brownian_particles):
            x_bm_track = self.data.T[self.brownian_particles_index[brownian_particle_i]][0][[i for i in range(frame_i + 1)]]
            y_bm_track = self.data.T[self.brownian_particles_index[brownian_particle_i]][1][[i for i in range(frame_i + 1)]]
            x_bm_track, y_bm_track = mask_periodic_plot(x_bm_track, y_bm_track)
    
            self._brownian_motion_track[brownian_particle_i].set_data(x_bm_track, y_bm_track)
    
    return self._gas_state, self._brownian_particles_state, self._brownian_motion_track


def animate_vel_distribution(
        self,
        frame_i: int
        ) -> Figure:
    """
    Generates the frames for the velocity distribution 
    animation.

    Returns a figure of velocity historgram.

    
    Keyword arguments:
        
    frame_i (int) -- the index of the i-th frame 
    """
    
        
    vel_cum_count = np.sum(self.vel_count[0:frame_i + 1], axis = 0)/(frame_i + 1)
    for bar, height in zip(self._vel_hist, vel_cum_count[::10]):
        bar.set_height(height)
        
    return self._vel_hist
    

def gas_animation_pbc(
        self, 
        exit_file_name: str, 
        fps: int, 
        dpi: int
        ) -> None:
    
    """
    Generates the animations of the gas with periodic boundary 
    conditions and saves as a .mp4 file.

    Keyword arguments:
        
    exit_file_name (str) -- exit file name
    fps (int) -- frames per second
    dpi (float) -- dots per inch 
    
    """

    fig = plt.figure()
    ax = plt.axes(xlim = (0., self.box_width), ylim = (0., self.box_width))
    self._gas_state, = ax.plot([], [], 
                              marker = "o", 
                              linestyle = " ", 
                              markerfacecolor = "cornflowerblue", 
                              markeredgecolor = "cornflowerblue", 
                              markersize=5
                              )

    self._brownian_particles_state, = ax.plot([], [], 
                                            marker = "o", 
                                            linestyle = " ", 
                                            markerfacecolor = "orange", 
                                            markeredgecolor = "orange", 
                                            markersize=9
                                            )

    self._brownian_motion_track = [ax.plot([], [], color = "tan")[0] for _ in range(self.n_brownian_particles)]

    molecular_dynamics_animation = animation.FuncAnimation(
        fig = fig, 
        func = functools.partial(animate_gas,self), 
        init_func = functools.partial(init_state, self), 
        frames = [i for i in range(self.total_time_steps)][::10]
        )
    plt.close()

    molecular_dynamics_animation.save(
        exit_file_name,
        fps = fps, 
        dpi = dpi
        )
    

def gas_animation_bhw(
        self,
        exit_file_name: str, 
        fps: int, 
        dpi: int
        ) -> None:
    
    """
    Generates the animations of the gas with bottom 
    hard wall boundary conditions and saves as a .mp4 file.

    Keyword arguments:
        
    exit_file_name (str) -- exit file name
    fps (int) -- frames per second
    dpi (float) -- dots per inch 
    
    """

    x_wall = np.array([i*self.box_width/10. for i in range(11)])

    fig = plt.figure()
    ax = plt.axes(xlim = (0., self.box_width), ylim = (-2., 1.5*self.box_width))

    plt.plot(x_wall, 0.*x_wall, "k")
    for i in range(-1,41):
        x_aux = np.array([j*self.box_width/40. + i*self.box_width/40. for j in range(3)])
        y_aux = 2.*x_aux - 2.*(0.5*i + 1)
        plt.plot(x_aux, y_aux, ":k")

    self._gas_state, = ax.plot([], [], 
                              marker = "o", 
                              linestyle = " ", 
                              markerfacecolor = "cornflowerblue", 
                              markeredgecolor = "cornflowerblue", 
                              markersize=5
                              )

    self._brownian_particles_state, = ax.plot([], [], 
                                            marker = "o", 
                                            linestyle = " ", 
                                            markerfacecolor = "orange", 
                                            markeredgecolor = "orange", 
                                            markersize=9
                                            )

    self._brownian_motion_track = [ax.plot([], [], color = "tan")[0] for _ in range(self.n_brownian_particles)]

    molecular_dynamics_animation = animation.FuncAnimation(
        fig = fig, 
        func = functools.partial(animate_gas,self), 
        init_func = functools.partial(init_state, self), 
        frames = [i for i in range(self.total_time_steps)][::10]
        )
    plt.close()

    molecular_dynamics_animation.save(
        exit_file_name,
        fps = fps, 
        dpi = dpi
        )
    

def vel_distribution_animation(
        self, 
        exit_file_name: str
        ) -> None:
    
    """
    Generates the animations of the velocity distribution 
    as an histogram and saves as a .mp4 file.

    Keyword arguments:
        
    exit_file_name (str) -- exit file name

    """

    fig, ax= plt.subplots()
    ax.set(xlim = (0, self.vel_bins.max()), ylim = (0, 1.05*(self.vel_count[self.total_time_steps - 1].max())))
    plt.title("Distribution of the molecules's velocity modulus")
    plt.xlabel("velocity")
    plt.ylabel("density")
    self._vel_hist = ax.bar(
        x = self.vel_bins[::10], 
        height = np.empty(len(self.vel_bins[::10])), 
        width = 0.085,
        color = "cornflowerblue"
        )

    vel_dist_animation = animation.FuncAnimation(
        fig = fig, 
        func = functools.partial(animate_vel_distribution, self), 
        frames = [i for i in range(self.total_time_steps)][::10]
        )
        
    plt.close()

    vel_dist_animation.save(
        exit_file_name,
        fps=50, 
        dpi = 150
        )
