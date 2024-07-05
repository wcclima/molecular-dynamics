# Molecular Dynamics
Repo for the molecular dynamics project.

## 1 - Objective

The goal of this project is to simulate a gas by modelling its molecules as classical particles interacting with each other via the Lennard-Jones potential. The simulation allows us to produce animations depicting the gas molecules's movement, to see the gas achieve its thermal equilibrium, to measure its temperature and to change its temperature to show phase transition. 

## 2 - Repo organisation

**`notebooks`: Notebooks demonstrating the simulations**
- `molecular_dynamics_notebook.ipynb` : Notebook discussing different aspects of the simulation

**`modules`: Simulation architecture**
- `molecular_dynamics.py` :

## 3 - The model

### 3.1 - Interaction potential

Here we will consider a set of $N$ *classical* molecules of mass $m$ that interact pair-wise via the Lennard-Jones potential and moving according Newton's laws in two dimensions. The Lennard-Jones potential is given by 

$$
V(r_{ij}) = 4\epsilon\left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6\right],
$$

where $r_{ij} = r_{ji}$ is the Euclidean distance between the molecules $i$ and $j$, $\epsilon$ is the strength of the potential and $\sigma$ is a length scale.

<p align="center">
  <img src="https://github.com/wcclima/molecular-dynamics/blob/main/lennard_jones.png" />
</p>

We see from the shape of the Lennard-Jones potential above that for long distances, i.e. for $r \gtrsim \sigma$, the potential is attractive, while for $r \lesssim \sigma$ the potential is strongly repulsive. The potential also has a minimum at $r_0 = 2^\frac{1}{6}\sigma$, where $V(r_0) = -\epsilon$.

This potential is a good model for innert gases, like ideal gases, that cannot for chemical bonds, i.e. cannot exchange electrons. The $\left(\frac{\sigma}{r}\right)^6$ tail models dipole interactions that are produce by the interaction of the molecules electronic clouds via the electorstatic force, that distorts the clouds, attractive producing electric dipole forces. The $\left(\frac{\sigma}{r}\right)^{12}$-term models the strong repulsion that appears out of the superpposition of the molecules electronic clouds, and it is produced by Pauli's exclusion principle. This principle express the empirical fact that two electrons cannot occupy the same state.


### 3.2 - Boundary conditions

**Periodic boundary conditions:** For simplicit, we will consider periodic boundary conditions. This means that in a lattice of size $L$, if the molecule moves to a position $x = L + a$, the molecule is the moved to the position $x' = x - L = a$, assuming here $a < L$. Hence, if the molecules leaves the lattice on the right, it reappers on the left. And vice-versa. If the molecule moves to the position $x = a < 0$, then it is moved to the position $x' = a + L < L$.

This choice of boundary conditions have some advantages. In a real gas we typically have $N \sim 10^{23}$ molecules, which would be very costly to simmulate in a computer. For small set of molecules, say $N \sim 100$, assuming the the boundary of the lattice are hard walls as in a vessel might produce undesirable boundary effects that one would not find in real-life situations. The periodic boundary condition then allows us to simulate the situation in which we are observing a sample of the bulk of the gas, away from the vessel walls at a constant density.

**Bottom hard wall condition:**

## 4 - Method

### 4.1 - Numerical method

To numerically compute the position of the gas molecules, we will use the so-called **Verlet method**, see e.g. Ref. [1]. This method is suitable when the interaction potential does not depend on the velocity. Moreover it is a method easy to implement and has good numerical precision.

Let us consider a function $y(t)$ satisfy the equation of motion $\frac{d^2y(t)}{dt^2} = f[y(t)]$ and its Taylor expansion going forward and backward in time around the instant $t_i$ with a time step $\Delta t$,

$$
y(t_i + \Delta t) = y(t_i) + \Delta t\frac{dy(t_i)}{dt} + \frac{1}{2}\Delta t^2\frac{d^2y(t_i)}{dt^2} + \frac{1}{6}\Delta t^3\frac{d^3y(t_i)}{dt^3} + O(\Delta t^4)
$$
and
$$
y(t_i - \Delta t) = y(t_i) - \Delta t\frac{dy(t_i)}{dt} + \frac{1}{2}\Delta t^2\frac{d^2y(t_i)}{dt^2} - \frac{1}{6}\Delta t^3\frac{d^3y(t_i)}{dt^3} + O(\Delta t^4).
$$

Denoting $y_i \equiv y(t_i)$, $y_{i-1} \equiv y(t_i - \Delta t)$, $y_{i+1} \equiv y(t_i + \Delta t)$ and $f(y) = $, we have that

$$
y_{i + 1} = 2 y_i - y_{i - 1} + \Delta t^2f(y_i) + O(\Delta t^4),
$$

hence the error is of order $O(\Delta t^4)$. The Verlet method thus approximates $\frac{d^2y}{dt^2}$ symmetrically using differences centered at the instant $t_i$ with high accurace. We can also obtain the velocity $v_i = \frac{dy_i}{dt}$ from the Verlet method by taking

$$
y_{i+1} - y_{i-1} = 2\Delta t v_i + O(\Delta t^3)\hspace{1cm} \rightarrow \hspace{1cm} v_i = \frac{y_{i+1} - y_{i-1}}{2\Delta t} + O(\Delta t^3).
$$

### 4.2 - Initial conditions

For the initial condition, e assume that the gas has a fixed density $\rho$, with its molecules starting on the vertices of a square lattice of side $L$ and cell size $d = L/\sqrt{\rho}$. We give the option of adding a random displacement $s = \alpha d$ both in the $x$- and $y$-directions, with $\alpha$ random variable homogeneously distributed in the interval $[-0.25, 0.25]$. For the velocity $\vec{v}$, we assume that $v_x$ and $v_y$ components are Gaussian random variables.

<p align="center">
  <img src="https://github.com/wcclima/molecular-dynamics/blob/main/example_initial_conditions.png" />
  <em>Example of initial condition in a lattice of side L = 20, with the arrows displaying the initial velocities of the molecules. The density was assumed 0.25.</em>
</p>

Given the initial positions and velocities of the molecules, we compute the positions at the previous step, needed to start the Verlet method, by simply taking the positions for the $i$-th molecule at the step as

$$
y_i(t_{-1}) = y_i(t_0) - v_i(t_0)\Delta t.
$$

### 4.3 - Computation of the forces

The expression for the interaction force between two molecules is

$$
\vec{F}(\vec{x} - \vec{x}') = \frac{24 \epsilon}{\sigma}\left[2\left(\frac{\sigma}{r}\right)^{13} - \left(\frac{\sigma}{r}\right)^7\right]
\frac{\vec{x} - \vec{x}'}{r} 
$$

with $r = \| \vec{x} - \vec{x}'\|$. 

To simplify the numerical computation we will only consider the interaction between molecules that are apart by up to a distance of $3\sigma$. This introduces an error in the computation of force of the order $1/3^7 \sim 10^{-4}$. Hence, if we choose $\Delta t \sim 10^{-3}$, the error in the equation for the position is still consistent with the error in the equation for the velocity.

We also need to take into account the periodic boundary condition when computing the force so the reappearnce of a molecule on the other side of the lattice is consistent with the dynamics. This is accomplished by computing the shortest distance between two molecules and using it to compute the force between them. Hence, suppose two molecules $i$ and $j$ at positions $x_i$ and $x_j$. If their coordinate distance is less than $L/2$, then we simply use $x_i - x_j$, otherwise we have to use the complement of the coordinate distance given by $x_i - x_j - \textrm{sgn}(x_i - x_j)L$.


<p align="center">
  <img src="https://github.com/wcclima/molecular-dynamics/blob/main/example_initial_forces.png" />
  <em>Example of initial condition in a lattice of side L = 20, with the arrows now displaying the initial total forces on the molecules. The density was assumed 0.25.</em>
</p>

## 5 - Statistical measurements

### 5.1 - Maxwell-Boltzmann distribution

Here we will consider an ideal or diluted gas, i.e. a gas whose the interaction energy bewteen the molecules is negligible when compared to the kinetic energy. This is done by choosing the molecules's initial state such that their separation is $\sim \sigma$ and velocity sufficiently high. Note that the initial state will not correspond to a thermal equilibrium state, so the collisions between the molecules are still important to take the system to equilibrium.

Hence, let us assume that our system is well approximated by an ideal gas. Since the probability of finding the gas in a state with energy between $E$ and $E + dE$ is given $p(E)dE \propto \textrm{e}^{-\frac{E}{kT}}dE$, where $k$ is Boltzmann's constant and $T$ the temperature. For an ideal gas the particles are independent from each other and the energy is the kinetic energy, $E = \frac{1}{2}mv^2$. Hence, the propability of finding an atom with velocity modulus between $v$ and $v + dv$ in 2 spatial dimentions is  

$$
\eqalign{
p(E)dE & \propto \exp{\left(-\frac{m}{2kT}v^2\right)}\vec{v}\cdot d\vec{v} \cr 
& \propto v\exp{\left(-\frac{m}{2kT}v^2\right)}dv.
}
$$

The normalised propbability distribution for the velocity modulus then is

$$
f(v)dv = \frac{m}{kT} v \exp{\left(-\frac{m}{2kT}v^2\right)}dv,
$$

the Maxwell-Boltzmann distribution.

### 5.2 - Measuring the temperature

To measure the temperature of our gas once it reaches thermal equilibrium, we can measure the velocity modulus distribution and then fit the resulting curve. For the fitting, we notice that for the discrete probability distribution $g(v)$ we have

$$
\ln g(v) = \ln\left(\frac{m\Delta v}{kT}\right) + \ln v - \frac{m}{2kT}v^2,
$$

and thus

$$
\ln\left(\frac{g(v)}{v}\right) = \ln\left(\frac{m\Delta v}{kT}\right) - \frac{m}{2kT}v^2,
$$

where $\Delta v$ is the finite size of the velocity bins. Thus, we can simply use a linear regression, after a log-log transformation of the velocity distribution.

For a small system such as the ones we are able to simulate, there will be a considerable amount of statistical fluctuation in the velocity modulus populations. A way to conter-balance these fluctuations is to average the velocity distribution over time.

Another possibility is to use the equipartion theorem, where each degree of freedom corresponds to $\frac{1}{2}kT$ in the mean energy. For a free particle in 2 dimension, we have that

$$
\left\langle E \right\rangle = kT
$$

and

$$
kT = \frac{m}{2}\left\langle v^2 \right\rangle. 
$$

Thus an alternative way to measure the temperature is to measure the mean of the squared velocity.

## 6 - Thermostat

## 7 - References

1. N. J. Giordano and H Nakanishi, *Computational Physics* (Pearson Prentice Hall, New Jersey, 2006).
2. R.H. Landau, M.J. Páez and C.C. Bordeianu, *Computation Physics: Problem Solving with Python* (Wiley-VCH, Weinheim, 2015). 
