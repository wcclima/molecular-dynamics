
import numpy as np

import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

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
    """

    self.vel_bins = np.array([i*bin_size + 0.5*bin_size for i in range(bins)])
    self.vel_count = np.empty((self.total_time_steps, bins))
    v_count = np.zeros(bins)

    for frame_i in range(self.total_time_steps):
        vx, vy = self.data[frame_i][2:4]

        v_count = velocity_modulus_fractional_count(
            vx[self.gas_molecules_index], 
            vy[self.gas_molecules_index], 
            self.n_gas_molecules, 
            bins, 
            bin_size
        )

        self.vel_count[frame_i] = v_count


def measure_kT(
        self, 
        method: str
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

    Returns the value value of kT (float).

    Keyword arguments:
    
    method (str) -- the kT measurement method:
    'instantaneous fitting', 'cumulative fitting' 
    and 'equipartition'.

    """

        
    if method == "instantaneous fitting":
        remove_zeros_mask = self.vel_count[-1] > 0.

        X = (self.vel_bins[remove_zeros_mask]**2).reshape(-1,1)
        y = np.log((self.vel_count[-1]/self.vel_bins)[remove_zeros_mask])

        lin_reg = LinearRegression()
        lin_reg.fit(X,y)
        lin_reg.coef_, lin_reg.intercept_ 

        dv = 2.*self.vel_bins[0]

        self.average_kT = 0.5*(np.exp(-lin_reg.intercept_)*dv - 1/(2.*lin_reg.coef_[0]))

        print(f"Measured kT after {self.total_time_steps} steps from fitting: {self.average_kT}")
        print(" ")

        fitted_vel_count = -2.*lin_reg.coef_[0]*self.vel_bins*np.exp(lin_reg.coef_[0]*self.vel_bins**2)*dv

        fig, ax = plt.subplots(dpi = 150)
        plt.plot(self.vel_bins, self.vel_count[-1])
        plt.plot(self.vel_bins, fitted_vel_count)
        plt.show()

    elif method == "cumulative fitting":

        cum_vel_count = np.mean(self.vel_count, axis = 0)
        remove_zeros_mask = cum_vel_count > 0.

        X = (self.vel_bins[remove_zeros_mask]**2).reshape(-1,1)
        y = np.log((cum_vel_count/self.vel_bins)[remove_zeros_mask])

        lin_reg = LinearRegression()
        lin_reg.fit(X,y)
        lin_reg.coef_, lin_reg.intercept_ 

        dv = 2.*self.vel_bins[0]

        self.average_kT = 0.5*(np.exp(-lin_reg.intercept_)*dv - 1/(2.*lin_reg.coef_[0]))

        print(f"Measured kT after {self.total_time_steps} steps from fitting: {self.average_kT}")
        print(" ")

        fitted_vel_count = -2.*lin_reg.coef_[0]*self.vel_bins*np.exp(lin_reg.coef_[0]*self.vel_bins**2)*dv

        fig, ax = plt.subplots(dpi = 150)
        plt.plot(self.vel_bins, cum_vel_count)
        plt.plot(self.vel_bins, fitted_vel_count)
        plt.show()

    elif method == "equipartition":

        vx, vy = self.data[-1][2:4]
        v_squared_mean = np.mean(vx*vx + vy*vy)

        print(f"Measured kT after {self.total_time_steps} steps from the equipartition theorem: {0.5*v_squared_mean}")

    else:
        raise WrongTemperatureMeasurementMethod(
            f"Error: '{method}' is not a method."
        )
