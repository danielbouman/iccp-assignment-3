"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def animate_wavefunction(probability,L,a,duration):
    fig, ax = plt.subplots()

    x = np.arange(0, L, a)        # x-array
    line, = ax.plot(x, np.sin(x))

    def animate(i):
        line.set_ydata(probability[:,i])  # update the data
        return line,

    #Init only required for blitting to give a clean slate.
    def init():
        line.set_ydata(np.ma.array(x, mask=True))
        return line,

    ani = animation.FuncAnimation(fig, animate, np.arange(1, duration), init_func=init,
        interval=35, blit=True,repeat=False)
    plt.axis([0,L,-0.05,0.5])
    potential = np.zeros((L/a,1),dtype=float)
    potential[100] = 0.3
    plt.plot(x,potential)
    plt.show()
    return