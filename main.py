# ICCP assignment 3
from QD.dynamics import CrankNicolson
import numpy as np

# Parameters
a = 0.5         # Spatial resolution
L = 50          # Square box size
sigma_x = 2.5   # Wavefunction shape x-direction
sigma_y = 2.5   # Wavefunction shape y-direction
mu_x = 25       # Wavefunction position x-direction
mu_y = 5        # Wavefunction position y-direction
k_x = 1.1*(2*np.pi)/a         # Wave vector x-direction
k_y = 0*(2*np.pi)/a       # Wave vector y-direction

# Time evolution
tau = 0.1
duration = 200

particle = CrankNicolson(a,L,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y) # Initialize particle and momentum operators
# particle.potential("double slit",25,100,1,5)                    # Initialize potential
particle.potential("harmonic trap",0.2)                    # Initialize potential
# particle.timeEvolution(tau,duration)                            # Start time evolution of the particle
# particle.saveData("")
# particle.plot("animate",saveAnimation=True)                     # Plot the result