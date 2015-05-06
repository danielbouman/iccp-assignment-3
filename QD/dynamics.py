import numpy as np
from numpy import linalg as LA
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import sys
from anim2D import animate_wavefunction
# import time
# import pylab
# from mpl_toolkits.mplot3d import Axes3D

class CrankNicolson:
  
  def __init__(self,a,L,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y):
    self.gridLength = int(L/a) # box length
    self.a = a # spatial resolution
    self.grid1D = np.linspace(0,L,self.gridLength)
    self.L = L

    # Define the momentum part of the Hamiltonian matrices
    c = 1/(a**2)
    xNeighbors = sp.diags([c,-4*c,c],[-1,0,1],shape=(self.gridLength,self.gridLength))
    self.H = sp.kron(sp.eye(self.gridLength),xNeighbors) + sp.diags([c,c],[-self.gridLength,self.gridLength],shape=(self.gridLength**2,self.gridLength**2))
    
    # Wave function
    psi_x = np.multiply(np.exp((-1/sigma_x)*np.power(self.grid1D-mu_x,2)),np.exp(-1j*k_x*self.grid1D))
    psi_y = np.multiply(np.exp((-1/sigma_y)*np.power(self.grid1D-mu_y,2)),np.exp(-1j*k_y*self.grid1D))
    self.psi = np.outer(psi_y,psi_x).flatten()
    
  def potential(self,function,*args):
    V = sp.lil_matrix((self.gridLength,self.gridLength))
    if str.lower(function) == "wall":
      """
      Horizontal wall
      arg 0 : y-position of wall
      arg 1 : height of wall
      """
      V[args[0]/self.a,:] = args[1]

    if str.lower(function) == "double slit":
      """
      Double slit wall
      arg 0 : y-position of wall
      arg 1 : height of wall
      arg 2 : slit distances from the center
      arg 3 : width of slits
      """
      V[args[0]/self.a,:] = args[1]
      endIndexSlitLeft = int(self.gridLength/2 - args[2]/self.a)
      startIndexSlitRight = int(self.gridLength/2 + args[2]/self.a)
      V[args[0]/self.a,int(endIndexSlitLeft-args[3]/self.a):endIndexSlitLeft] = 0
      V[args[0]/self.a,startIndexSlitRight:int(startIndexSlitRight+args[3]/self.a)] = 0
    
    self.H = self.H + sp.diags(V.reshape((1,self.gridLength**2)).toarray(),[0])
      

  def normalize_wavefunction(self):
    self.psi = (1/np.linalg.norm(self.psi))*self.psi
    
  def timeEvolution(self,tau,duration):
    self.duration = duration
    # Define A and B matrices
    A = sp.identity(self.gridLength**2) - tau/(1j)*self.H
    B = sp.identity(self.gridLength**2) + tau/(1j)*self.H
    
    self.time_evolved_psi = np.zeros((self.gridLength**2,duration),dtype=complex)

    # Start time evolution of particle
    for i in range(0,duration):
      # print(i)
      # Solve linear equation A*psi(t + tau) = B*psi(t)
      self.time_evolved_psi[:,i],_ = linalg.bicgstab(A,B.dot(self.psi).transpose())
      self.psi = self.time_evolved_psi[:,i]

  def plot2D(self,plotStyle):
    time_evolved_probability = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi))).reshape(self.gridLength,self.gridLength,self.duration)
    x,y = np.meshgrid(self.grid1D,self.grid1D)

    if plotStyle == "plot":
        for i in range(0,self.duration):
          grid = time_evolved_probability[:,:,5*i]
          plt.imshow(grid, interpolation='none')
          plt.show()
    else:
        fig = plt.figure()
        ax = plt.axes(xlim=(0, self.L), ylim=(0, self.L))  

        # animation function
        def animate(i): 
            z = time_evolved_probability[:,:,i]
            cont = plt.contourf(x, y, z,9)
            return cont  

        anim = animation.FuncAnimation(fig, animate, interval= 200,  repeat_delay=1000, frames=self.duration)
        plt.show()


