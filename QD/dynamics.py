import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from anim import animate_wavefunction


class CrankNicolson:
  
  def __init__(self,a,L,waveform="gaussian",*args):
    self.waveform = waveform
    self.gridLength = int(L/a +1)
    self.a = a # spatial resolution
    self.grid = np.linspace(0,L,self.gridLength)
    
    # Define the momentum operator matrix
    c = -1/(a**2)
    self.H = sp.diags([c,-2*c,c],[-1,0,1],shape=(self.gridLength,self.gridLength))
    
    # Wave function
    if waveform == "gaussian":
      """
      arg 0 : sigma
      arg 1 : k
      arg 2 : mu
      """
      self.psi = np.multiply(np.exp((-1/args[0]**2)*np.power(self.grid-args[2],2)),np.exp(-1j*args[1]*self.grid))
    if waveform == "eigenstate":
      """
      arg 0 : n
      """
      self.psi = np.sin(args[0]*np.pi/L*self.grid)
    # Normalize wave function
    self.psi = np.multiply((1/(np.linalg.norm(self.psi))),self.psi)

    
  def potential(self,function="infinite square well",*args):
    self.potentialName = function
    # Initialize potential array
    self.V = np.zeros(self.gridLength)
    if str.lower(function) == "rectangular barrier":
      """
      Rectangular barrier
      arg 0 : position of barrier
      arg 1 : height of barrier
      arg 2 : width of barrier
      """
      self.barrierStart = args[0]/self.a
      self.barrierEnd = args[0]/self.a+args[2]/self.a
      # Add barrier to potential array
      self.V[self.barrierStart:self.barrierEnd] = args[1]
      # Add potential to Hamiltonian operator
      self.H = self.H + sp.diags([self.V],[0])
    
  def timeEvolution(self,tau,duration):
    self.duration = duration
    # Define A and B matrices
    A = sp.identity(self.gridLength) - tau/(2j)*self.H
    B = sp.identity(self.gridLength) + tau/(2j)*self.H
    self.time_evolved_psi = np.zeros((self.gridLength,duration),dtype=complex)
    # Start time evolution
    for i in range(0,duration):
      self.time_evolved_psi[:,i],_ = linalg.bicgstab(A,B.dot(self.psi).transpose(),tol=1e-10)
      self.psi = self.time_evolved_psi[:,i]
    if str.lower(self.potentialName) == "rectangular barrier":
      return np.trapz(np.abs(self.psi[self.barrierEnd:])**2)
    
  def animate(self,saveAnimation=False):
    # Time evolved probability
    plotProb = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi)))
    fig, ax = plt.subplots()
    if saveAnimation == True:
      # Set up formatting for the movie files
      Writer = animation.writers['ffmpeg']
      writer = Writer(fps=15, metadata=dict(artist='Bouman and Goodenough'))

    line, = ax.plot(self.grid, np.sin(self.grid))

    def animate(i):
        line.set_ydata(plotProb[:,i])  # update the data
        return line,

    #Init only required for blitting to give a clean slate.
    def init():
        line.set_ydata(np.ma.array(self.grid, mask=True))
        return line,

    ani = animation.FuncAnimation(fig, animate, np.arange(1, self.duration), init_func=init,
        interval=35, blit=True,repeat=False)
    plt.axis([0,self.gridLength*self.a,0,0.15])
    plt.plot(self.grid,self.V*0.01)
    if saveAnimation == True:
      print("Saving animation...")
      start = time.time()
      ani.save(self.potentialName+'.mp4', writer=writer)
      end = time.time()
      print(end - start)
      print("Done. Animation saved as "+self.potentialName+".mp4")
    else:
      plt.show()
    return