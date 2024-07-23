
import numpy as np

import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

from typing import Optional

class WrongTemperatureMeasurementMethod(Exception):
    def __init__(self, message):
        super().__init__(message)


def velocity_modulus_fractional_count(
        vx: np.ndarray, 
        vy: np.ndarray, 
        samples: int, 
        bins: int, 
        bin_size: float#
        ) -> np.ndarray:
    """
    Counts the fractional number of molecules with 
    velocity modulus within a certain velocity bin dv.

    Returns an array with the fraction of molecules with 
    velocity modulus within v and v + dv for a certain 
    velocity range.

    
    Keyword arguments:
    
    vx (np.ndarray) -- the velocity x component
    vy (np.ndarray) -- the velocity y component
    samples (int) -- number of velocity samples
    bins (int) -- number of histogram bins
    bin_size (float) -- size of the histogram bins
    """
    
    v = np.sqrt(vx*vx + vy*vy)
    count_list = []
    for i in range(bins):
        count_list.append(np.sum(np.logical_and(v > i*bin_size, v < (i + 1)*bin_size)))
    
    return np.array(count_list)/samples


def vel_squared_mean(
        vx: np.ndarray, 
        vy: np.ndarray
        ) -> float:
    """
    Computes the instantaneous mean squared 
    velocity.
    
    Returns the mean squared velocity value.

    Keyword arguments:
    
    vx (np.ndarray) -- the velocity x component
    vy (np.ndarray) -- the velocity y component
    """
    
    v_squared = vx*vx + vy*vy
    return np.mean(v_squared)


def vel_distribution(
        self, 
        bins: int, 
        bin_size: float
        ) -> None:
    """
    Generates the fractional count of the velocity modulus
    for each time step of the gas history.

    
    Keyword arguments:
    
    bins (int) -- number of histogram bins
    bin_size (float) -- size of the histogram bins
    smoothing_window (int) -- size of the smoothing window, in steps
    """

    self._vel_bins = np.array([i*bin_size + 0.5*bin_size for i in range(bins)])
    self._vel_frac_count = np.empty((self.total_time_steps, bins))
    v_frac_count = np.zeros(bins)

    for frame_i in range(self.total_time_steps):
        vx, vy = self.data[frame_i][2:4]

        v_frac_count = velocity_modulus_fractional_count(
            vx[self.gas_molecules_index], 
            vy[self.gas_molecules_index], 
            self.n_gas_molecules, 
            bins, 
            bin_size
        )

        self._vel_frac_count[frame_i] = v_frac_count


def time_average_func(
        a: np.ndarray,
        average_window: int,
        frame_i: int
        ):

    window_size = np.min([average_window, frame_i + 1])

    a_averaged = np.sum(a[frame_i + 1 - window_size:frame_i + 1], axis = 0)
    a_averaged = a_averaged/window_size


    return a_averaged


def measure_kT(
        self, 
        method: str,
        average_window: Optional[int]
        ) -> None:
    """
    Computes the value of k times the temperature,
    where k is Boltzmann's constant by fitting the 
    velocity distribution using the instantaneous 
    Maxwell-Boltzmann distribution, by fitting the 
    velocity distribution using the time-averaged 
    Maxwell-Boltzmann distribution or by using 
    computing the instantaneous average kinectic 
    energy ans using the equipartition theorem. 


    Keyword arguments:
        method (str, values = ['instantaneous fitting', 'averaged fitting', 'equipartition']): 
            The kT measurement method.

        average_window (int, optional):
            The size of the time average window when method is 'averaged fitting'. 

    """

        
    if method == "instantaneous fitting":
        remove_zeros_mask = self._vel_frac_count[-1] > 0.

        X = (self._vel_bins[remove_zeros_mask]**2).reshape(-1,1)
        y = np.log((self._vel_frac_count[-1]/self._vel_bins)[remove_zeros_mask])

        lin_reg = LinearRegression()
        lin_reg.fit(X,y)
        lin_reg.coef_, lin_reg.intercept_ 

        dv = 2.*self._vel_bins[0]

        self.average_kT = 0.5*(np.exp(-lin_reg.intercept_)*dv - 1/(2.*lin_reg.coef_[0]))

        print(f"Measured kT after {self.total_time_steps} steps from fitting: {self.average_kT}")
        print(" ")

        fitted_vel_frac_count = -2.*lin_reg.coef_[0]*self._vel_bins*np.exp(lin_reg.coef_[0]*self._vel_bins**2)*dv

        fig, ax = plt.subplots(dpi = 150)
        plt.plot(self._vel_bins, self._vel_frac_count[-1])
        plt.plot(self._vel_bins, fitted_vel_frac_count)
        plt.show()

    elif method == "averaged fitting":

        vel_dist_averaged = time_average_func(self._vel_frac_count, average_window, self.total_time_steps)
        remove_zeros_mask = vel_dist_averaged > 0.

        X = (self._vel_bins[remove_zeros_mask]**2).reshape(-1,1)
        y = np.log((vel_dist_averaged/self._vel_bins)[remove_zeros_mask])

        lin_reg = LinearRegression()
        lin_reg.fit(X,y)
        lin_reg.coef_, lin_reg.intercept_ 

        dv = 2.*self._vel_bins[0]

        self.average_kT = 0.5*(np.exp(-lin_reg.intercept_)*dv - 1/(2.*lin_reg.coef_[0]))

        print(f"Measured kT after {self.total_time_steps} steps from fitting: {self.average_kT}")
        print(" ")

        vel_dist_fitted = -2.*lin_reg.coef_[0]*self._vel_bins*np.exp(lin_reg.coef_[0]*self._vel_bins**2)*dv

        fig, ax = plt.subplots(dpi = 150)
        plt.plot(self._vel_bins, vel_dist_averaged)
        plt.plot(self._vel_bins, vel_dist_fitted)
        plt.show()

    elif method == "equipartition":

        vx, vy = self.data[-1][2:4]
        v_squared_mean = np.mean(vx*vx + vy*vy)

        print(f"Measured kT after {self.total_time_steps} steps from the equipartition theorem: {0.5*v_squared_mean}")

    else:
        raise WrongTemperatureMeasurementMethod(
            f"Error: '{method}' is not a method."
        )
