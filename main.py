# ICCP assignment 3
from QD.dynamics import Particle

# Parameters
a = 0.5       # Spatial resolution
L = 100       # Domain size
sigma = 1     # Wavefunction shape
k = 5         # Wave vector
hbar = 1      # Reduced Planck constant
# Potential
pos = 50      # Potential position
amp = 10      # Potential height
<<<<<<< HEAD

# Time evolution
tau = 1
duration = 10


particle1 = Particle(a,L,sigma,k)     # Initialize particle
particle1.potential(pos,amp)          # Initialize potential
particle1.timeEvolution(tau,hbar,1)             # Start time evolution of particle

# particle1.plot()
