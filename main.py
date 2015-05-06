# ICCP assignment 3
from QD.dynamics import CrankNicolson

# Parameters
a = 0.5        # Spatial resolution
L = 40        # Square box size
sigma_x = 2.5   # Wavefunction shape x-direction
sigma_y = 2.5   # Wavefunction shape y-direction
mu_x = 20        # Wavefunction position x-direction
mu_y = 20       # Wavefunction position y-direction
k_x = 0     # Wave vector x-direction
k_y = 0     # Wave vector y-direction
# Potential
xStart =  20
xEnd = 25
yStart = 0
yEnd = 40
amp = 0         # Potential height

# Time evolution
tau = 0.1
duration = 100

particle = CrankNicolson(a,L,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y)     # Initialize particle
particle.normalize_wavefunction()      # normalize wavefunction so that probability sums up to unity
# particle.potential("wall")          # Initialize potential
particle.timeEvolution(tau,duration)             # Start time evolution of particle
particle.plot2D()
# particle1.plot()
