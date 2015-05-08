# ICCP assignment 3
import numpy as np
from QD.dynamics import CrankNicolson

# Parameters
a = 0.5                 # Spatial resolution
L = 100                 # Domain size
sigma = 9.5            # Wave function shape
k = 0.8*(2*np.pi)/a     # Wave vector
mu = 15                 # Wave function position

# Time evolution
tau = 0.1
duration = 140

pots = np.logspace(-2, 1, 100)
trans = np.zeros(len(pots))
i = 0

# particle = CrankNicolson(a,L,sigma,k,mu)     # Initialize particle
# particle.potential("rectangular barrier",40,50,2)          # Initialize potential
# _ = particle.timeEvolution(tau,duration)             # Start time evolution of particle
# particle.animate(saveAnimation=True)
for iii in pots:
  # print(iii)
  particle = CrankNicolson(a,L,sigma,k,mu)     # Initialize particle
  particle.potential("rectangular barrier",40,iii,2.6)          # Initialize potential
  trans[i] = particle.timeEvolution(tau,duration)             # Start time evolution of particle
  # print(iii)
  i = i + 1
print('Done.')