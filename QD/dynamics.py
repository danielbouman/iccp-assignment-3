import numpy as np
from numpy import linalg as LA
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import sys
# from anim2D import animate_wavefunction
# import time
# import pylab
# from mpl_toolkits.mplot3d import Axes3D

class CrankNicolson:
  
  def __init__(self,a,L,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y):
    self.gridLength = int(L/a) # box length
    self.a = a # spatial resolution
    self.grid1D = np.linspace(0,L,self.gridLength)
    
    # Define the momentum part of the Hamiltonian matrices
    c = 1/(a**2)
    xNeighbors = sp.diags([c,-4*c,c],[-1,0,1],shape=(self.gridLength,self.gridLength))
    self.H = sp.kron(sp.eye(self.gridLength),xNeighbors) + sp.diags([1,1],[-self.gridLength,self.gridLength],shape=(self.gridLength**2,self.gridLength**2))
    
    # Wave function
    psi_x = np.multiply(np.exp((-1/sigma_x)*np.power(self.grid1D-mu_x,2)),np.exp(-1j*k_x*self.grid1D))
    psi_y = np.multiply(np.exp((-1/sigma_y)*np.power(self.grid1D-mu_y,2)),np.exp(-1j*k_y*self.grid1D))
    self.psi = np.outer(psi_y,psi_x).flatten()
    
  def potential(self,function):
    if function == "wall":
      V = sp.csr_matrix((3,3),dtype="float")
      # V[] =
    
    # Add potential to hamiltonian matrix
    # xIndexStart = int(xStart/self.a)
    # xIndexEnd = int(xEnd/self.a)
    # yIndexStart = int(yStart/self.a)
    # yIndexEnd = int(yEnd/self.a)
    # self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] = self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] + amp

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

  def plot2D(self):
    time_evolved_probability = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi))).reshape(self.gridLength,self.gridLength,self.duration)
    # time_evolved_probability_y = np.real(np.multiply(self.time_evolved_psi_y,np.conj(self.time_evolved_psi_y)))
    # time_evolved_probabilit y= np.zeros((self.sizeX,self.sizeY,self.duration),dtype=float)

    # for i in range(0,self.duration):
        # time_evolved_probability[:,:,i] = np.outer(time_evolved_probability_x[:,i],time_evolved_probability_y[:,i])

    # print(time_evolved_probability[:,:,0])

    # animate_wavefunction(time_evolved_probability,self.L,self.a,self.duration)
    # plt.plot(self.xAxis,sel)

    x_mesh,y_mesh = np.meshgrid(self.grid1D,self.grid1D)
    
    # fig = plt.figure()
    # ax = plt.axes(xlim=(0, self.L_x), ylim=(0, self.L_y))
    # plt.xlabel(r'x')
    # plt.ylabel(r'y')

    # anim = animation.FuncAnimation(fig, animate, frames=Nt)

    # for i in range(0,self.duration):
    grid = time_evolved_probability[:,:,50]
    plt.imshow(grid, interpolation='none')
    plt.show()
