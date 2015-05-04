# Plot wavefunction through time in 1D

import numpy as np
import matplotlib.pyplot as plt
import time

def plot1D(time_evolved_psi,L,xAxis):
    time_evolved_probability = np.multiply(time_evolved_psi,np.conj(time_evolved_psi))

    plt.axis([0,L,0,6])
    plt.show()

    for ii in range(0,length(probability,1)):
        probability = time_evolved_probability[:,ii]
        plt.plot(xAxis,probability)
        plt.draw()
        time.sleep(0.5)
    return;