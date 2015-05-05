# ICCP assignment 3
from QD.dynamics import Particle

# Parameters
a = 0.25       # Spatial resolution
L_x = 40       # Domain size
L_y = 40 
sigma_x = 2.5     # Wavefunction shape
sigma_y = 6.5
k_x = 1/300         # Wave vector
k_y = -1/7         # Wave vector
mu_x = 10
mu_y = 30
hbar = 1      # Reduced Planck constant
# Potential
pos = 0      # Potential position
amp = 0      # Potential height

# Time evolution
tau = 0.1
duration = 500


particle1 = Particle(a,L_x,L_y,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y)     # Initialize particle
particle1.normalize_wavefunction()      # normalize wavefunction so that probability sums up to unity
particle1.potential(pos,amp)          # Initialize potential
particle1.timeEvolution(tau,hbar,duration)             # Start time evolution of particle
particle1.plot2D()
# particle1.plot()
