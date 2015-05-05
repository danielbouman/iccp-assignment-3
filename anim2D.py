"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def animate_wavefunction(probability,L,a,duration):
    

    # animation function
    def animate(i): 
        
        z = probability[i,:,0,:].T
        cont = plt.contourf(x, y, z, 25)
        if (tslice == 0):
            plt.title(r't = %1.2e' % t[i] )
        else:
            plt.title(r't = %i' % i)

        return cont  

    anim = animation.FuncAnimation(fig, animate, frames=Nt)

    anim.save('animation.mp4')

    return