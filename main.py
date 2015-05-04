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
duration = 10 # Amount of timesteps
=======
# Time evolution
tau = 1
>>>>>>> 7bde29dd2f6db2b06d90c111eefdef29be0c8882

particle1 = Particle(a,L,sigma,k)     # Initialize particle
particle1.potential(pos,amp)          # Initialize potential
particle1.timeEvolution(tau,hbar,1)             # Start time evolution of particle

# particle1.plot()
