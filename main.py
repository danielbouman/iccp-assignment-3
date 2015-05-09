# ICCP assignment 3
from QD.dynamics import CrankNicolson
import numpy as np

# Parameters
a = 0.5         # Spatial resolution
L = 40          # Square box size
sigma_x = 4.5   # Wavefunction shape x-direction
sigma_y = 4.5   # Wavefunction shape y-direction
mu_x = 20       # Wavefunction position x-direction
mu_y = 10        # Wavefunction position y-direction
k_x = -0.0*(2*np.pi)/a         # Wave vector x-direction
k_y = 1.00*(2*np.pi)/a       # Wave vector y-direction

# Time evolution
tau = 0.5
duration = 50

particle = CrankNicolson(a,L,tau,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y) # Initialize particle and momentum operators
# particle.potential("double slit",25,100,1,5)                    # Initialize potential
# particle.potential("harmonic",0.2)                    # Initialize potential
particle.timeEvolution(tau,duration)                            # Start time evolution of the particle
particle.wavefunctionComparison()
# particle.saveData("")
particle.plot("animate",saveAnimation=False)                     # Plot the result