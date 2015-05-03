# ICCP assignment 3
from QD.dynamics import Particle

# Parameters
a = 0.5       # Spatial resolution
L = 100       # Domain size
sigma = 1     # Wavefunction shape
k = 5         # Wave vector
# Dirac potential
pos = 50      # Potential position
amp = 10      # Potential height

particle1 = Particle(a,L,sigma,k)     # Initialize particle
particle1.potential(pos,amp)          # Initialize potential
particle1.timeEvolution()             # Start time evolution of particle

particle1