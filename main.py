# ICCP assignment 3
from QD.dynamics import Particle

# Parameters
a = 0.5       # Spatial resolution
L = 100       # Domain size
sigma = 9.5     # Wavefunction shape
k = 1         # Wave vector
mu = 10
hbar = 1      # Reduced Planck constant
# Potential
pos = 50      # Potential position
amp = 100      # Potential height

# Time evolution
tau = 0.1
duration = 500


particle1 = Particle(a,L,sigma,k,mu)     # Initialize particle
particle1.normalize_wavefunction()      # normalize wavefunction so that probability sums up to unity
particle1.potential(pos,amp)          # Initialize potential
particle1.timeEvolution(tau,hbar,duration)             # Start time evolution of particle
particle1.plot1D()
# particle1.plot()
