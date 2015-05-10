import numpy as np
from numpy import linalg as LA
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
  
class CrankNicolson:
  
  def __init__(self,a,L,wavefunction="gaussian",*args):
    self.gridLength = int(L/a + 1)
    self.a = a                                      # spatial resolution
    self.grid1D = np.linspace(0,L,self.gridLength)
    self.L = L                                      # Box length

    # Define the momentum part of the Hamiltonian matrices
    c = -1/(a**2)
    xNeighbors = sp.diags([c,-4*c,c],[-1,0,1],shape=(self.gridLength,self.gridLength))
    self.H = sp.kron(sp.eye(self.gridLength),xNeighbors) + sp.diags([c,c],[-self.gridLength,self.gridLength],shape=(self.gridLength**2,self.gridLength**2))
    # Create wave function
    if wavefunction == "gaussian":
      """
      Gaussian waveform
      args[0] : sigma_x
      args[1] : sigma_y
      args[2] : k_x
      args[3] : k_y
      args[4] : mu_x
      args[5] : mu_y
      """
      psi_x = np.multiply(np.exp((-1/args[0]**2)*np.power(self.grid1D-args[4],2)),np.exp(-1j*args[2]*self.grid1D))
      psi_y = np.multiply(np.exp((-1/args[1]**2)*np.power(self.grid1D-args[5],2)),np.exp(-1j*args[3]*self.grid1D))
    if wavefunction == "eigenstate":
      psi_x = np.sin(3*np.pi*self.grid1D/self.L)
      psi_y = np.sin(2*np.pi*self.grid1D/self.L)
    self.psi = np.outer(psi_y,psi_x).flatten()
    # Normalize wave function
    self.psi = (1/np.linalg.norm(self.psi))*self.psi
    

  def potential(self,function="infinite square well",*args):
    self.potentialName = function
    V = sp.lil_matrix((self.gridLength,self.gridLength))
    if str.lower(function) == "wall":
      """
      Horizontal wall
      arg 0 : y-position of wall
      arg 1 : height of wall
      """
      V[args[0]/self.a,:] = args[1]
      self.H = self.H + sp.diags(V.reshape((1,self.gridLength**2)).toarray(),[0])

    if str.lower(function) == "double slit" or str.lower(function) == "slit":
      """
      Double slit wall
      arg 0 : y-position of wall
      arg 1 : height of wall
      arg 2 : slit distances from the center
      arg 3 : width of slits
      """
      # Create hard wall potential
      V[args[0]/self.a,:] = args[1]
      endIndexSlitLeft = int(self.gridLength/2) - int(args[2]/self.a)
      startIndexSlitRight = int(self.gridLength/2) + int(args[2]/self.a)
      # Make slits in the wall
      V[args[0]/self.a,endIndexSlitLeft-int(args[3]/self.a)+1:endIndexSlitLeft+1] = 0
      V[args[0]/self.a,startIndexSlitRight:startIndexSlitRight+int(args[3]/self.a)] = 0
      # Reshape potential grid and add to Hamiltonian operator
      self.V = V
      self.H = self.H + sp.diags(V.reshape((1,self.gridLength**2)).toarray(),[0])
      
    if str.lower(function) == "harmonic trap" or str.lower(function) == "harmonic":
      """
      Harmonic trap potential
      arg 0 : omega
      """
      # Define x and y coordinates
      gridX = np.outer(np.ones((self.gridLength)),self.grid1D)
      gridY = np.outer(self.grid1D,np.ones((self.gridLength)))
      V = 0.5*(args[0]**2)*((gridX-gridX[(self.gridLength-1)/2,(self.gridLength-1)/2])**2+(gridY-gridY[(self.gridLength-1)/2,(self.gridLength-1)/2])**2)
      # Reshape potential grid and add to Hamiltonian operator
      self.H = self.H + sp.diags([V.flatten()],[0])
      self.potential = V
      self.potentialFlat = sp.diags([V.flatten()],[0])
    
  def timeEvolution(self,tau,duration):
    self.duration = duration
    # Define A and B matrices
    A = sp.identity(self.gridLength**2) - tau/(2j)*self.H
    B = sp.identity(self.gridLength**2) + tau/(2j)*self.H
    
    self.time_evolved_psi = np.zeros((self.gridLength**2,duration),dtype=complex)
    # Start time evolution of particle
    # Solve linear equation A*psi(t + tau) = B*psi(t)
    for i in range(0,duration):
      # print(i)
      self.time_evolved_psi[:,i],_ = linalg.bicgstab(A,B.dot(self.psi).transpose(),tol=1e-10)
      self.psi = self.time_evolved_psi[:,i]
    
  def wavefunctionComparison(self):
    probabilities = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi)))
    difference = np.zeros((self.gridLength**2,self.duration+1),dtype="complex")
    distance = np.zeros(self.duration+1)
    norm = np.zeros((self.duration+1),dtype="float")

    for i in range(self.duration):
        difference[:,i] = np.subtract(probabilities[:,i],probabilities[:,0])
        distance[i] = (self.a*2)*np.sum(np.real(np.multiply(difference[:,i],np.conj(difference[:,i]))))
        norm[i] = np.sum(probabilities[:,i])

    norm = np.subtract(norm,1)
    norm = np.multiply(norm,np.power(10,12))
    x = np.linspace(0,1,self.duration+1)
    plt.plot(x[0:-1],distance[0:-1])
    plt.show()
    plt.plot(x[0:-1],norm[0:-1])
    plt.show()

    # Save data
    np.savetxt('distance.dat', distance.flatten(), delimiter=',')


  def plot(self,plotStyle="",saveAnimation=False,fixColormap=False):
    time_evolved_probability = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi))).reshape(self.gridLength,self.gridLength,self.duration)
    x,y = np.meshgrid(self.grid1D,self.grid1D)

    if saveAnimation == True:
      # Set up formatting for the movie files
      Writer = animation.writers['ffmpeg']
      writer = Writer(fps=15, metadata=dict(artist='Bouman and Goodenough'))
    
    if str.lower(plotStyle) == "animate":
      print("Creating animation...")
      fig = plt.figure()
      ax = plt.axes(xlim=(0, self.L), ylim=(0, self.L))
      if fixColormap == True:
        maxProb = np.max(time_evolved_probability)
      else:
        maxProb = None
        
      def animate(i):
          cont = plt.contourf(x, y, time_evolved_probability[:,:,i],50,vmin=0,vmax=maxProb)
          return cont
          
      anim = animation.FuncAnimation(fig, animate, interval= 200,  repeat_delay=1000, frames=self.duration)
      
      if saveAnimation == True:
        print("Saving animation...")
        start = time.time()
        anim.save(self.potentialName+'.mp4', writer=writer)
        end = time.time()
        print(end - start)
        print("Done. Animation saved as "+self.potentialName+".mp4")
      else:
        plt.show()
        
    else:
      for i in range(0,self.duration):
        grid = time_evolved_probability[:,:,i]
        plt.imshow(grid, interpolation='none')
        plt.show()