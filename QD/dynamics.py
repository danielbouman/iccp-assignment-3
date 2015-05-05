import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from anim import animate_wavefunction


class Particle:
  
  def __init__(self,a,L,sigma,k,mu):
    self.L = L # box length
    self.a = a # spatial resolution
    self.xAxis = np.linspace(0,L,L/a)
    
    # Define the momentum operator matrix
    self.H = np.multiply(sp.eye(L/a,k=0),-2)
    self.H = self.H + sp.eye(L/a,k=-1)
    self.H = self.H + sp.eye(L/a,k=1)
    self.H = np.divide(self.H,-a**2).todense()
    
    # Wave function
    self.psi = np.multiply(np.exp((-1/sigma)*np.power(self.xAxis-mu,2)),np.exp(-1j*k*self.xAxis))

    
  def potential(self,pos,amp):
    # Add potential to hamiltonian matrix
    H_index = pos/self.a
    self.H[H_index,H_index] = self.H[H_index,H_index] + amp

  def normalize_wavefunction(self):
    inner_product = self.a*sum(np.multiply(self.psi,np.conj(self.psi)))
    self.psi = (1/np.sqrt(np.abs(inner_product)))*self.psi
    
  def timeEvolution(self,tau,hbar,duration):
    self.duration = duration
    # Define A and B matrices
    A = sp.identity(self.L/self.a) - tau/(1j*hbar)*self.H
    B = sp.identity(self.L/self.a) + tau/(1j*hbar)*self.H
    self.time_evolved_psi = np.zeros((self.L/self.a,duration),dtype=complex)
    # Time is run here
    for i in range(0,duration):
        self.time_evolved_psi[:,i],_ = linalg.bicgstab(A,B.dot(self.psi).transpose())
        self.psi = self.time_evolved_psi[:,i]
    
  def animate(self,interval=35):
    
    time_evolved_probability = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi)))
    fig, ax = plt.subplots()

    x = np.arange(0, self.L, self.a)        # x-array
    line, = ax.plot(x, np.sin(x))

    def animate(i):
        line.set_ydata(probability[:,i])  # update the data
        return line,

    #Init only required for blitting to give a clean slate.
    def init():
        line.set_ydata(np.ma.array(x, mask=True))
        return line,

    ani = animation.FuncAnimation(fig, animate, np.arange(1, self.duration), init_func=init,
        interval=interval, blit=True,repeat=False)
    plt.axis([0,self.L,-0.05,0.5])
    # potential = np.zeros((self.L/self.a,1),dtype=float)
    # potential[100] = 0.3
    # plt.plot(x,potential)
    plt.show()
    return