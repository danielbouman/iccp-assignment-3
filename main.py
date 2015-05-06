# ICCP assignment 3
from QD.dynamics import Particle

# Parameters
a = 0.25        # Spatial resolution
L_x = 40        # Domain length in x-direction
L_y = 40        # Domain lenght in y-direction
sigma_x = 2.5   # Wavefunction shape x-direction
sigma_y = 2.5   # Wavefunction shape y-direction
mu_x = 5        # Wavefunction position x-direction
mu_y = 20       # Wavefunction position y-direction
k_x = 99999     # Wave vector x-direction
k_y = 0         # Wave vector y-direction
# Potential
xStart =  20
xEnd = 21
yStart = 0
yEnd = 40
amp = 1         # Potential height

# Time evolution
tau = 0.1
duration = 100

particle1 = Particle(a,L_x,L_y,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y)     # Initialize particle
particle1.normalize_wavefunction()      # normalize wavefunction so that probability sums up to unity
particle1.potential(xStart,yEnd,yStart,yEnd,amp)          # Initialize potential
particle1.timeEvolution(tau,duration)             # Start time evolution of particle
particle1.plot2D()
# particle1.plot()
