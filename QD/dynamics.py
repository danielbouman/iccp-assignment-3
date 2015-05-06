import numpy as np
from numpy import linalg as LA
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
# from anim2D import animate_wavefunction
# import time
# import pylab
# from mpl_toolkits.mplot3d import Axes3D

class Particle:
  
  def __init__(self,a,L_x,L_y,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y):
    self.sizeX = int(np.ceil(L_x/a)) # box length
    self.sizeY = int(np.ceil(L_y/a)) # box length
    self.a = a # spatial resolution
    self.xAxis = np.linspace(0,L_x,L_x/a)
    self.yAxis = np.linspace(0,L_y,L_y/a)
    
    # Define the momentum part of the Hamiltonian matrices
    self.H_x = np.divide(sp.diags([1,-2,1],[-1,0,1], shape=(self.sizeX,self.sizeX)),-a**2).todense()
    self.H_y = np.divide(sp.diags([1,-2,1],[-1,0,1], shape=(self.sizeY,self.sizeY)),-a**2).todense()

    # Wave function
    psi_x = np.multiply(np.exp((-1/sigma_x)*np.power(self.xAxis-mu_x,2)),np.exp(-1j*k_x*self.xAxis))
    psi_y = np.multiply(np.exp((-1/sigma_y)*np.power(self.yAxis-mu_y,2)),np.exp(-1j*k_y*self.yAxis))
    self.psi = np.outer(psi_x,psi_y)
    
  def potential(self,xStart,xEnd,yStart,yEnd,amp):
    # Add potential to hamiltonian matrix
    xIndexStart = xStart/self.a
    xIndexEnd = xEnd/self.a
    yIndexStart = yStart/self.a
    yIndexEnd = yEnd/self.a
    self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] = self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] + amp
    self.H_y[yIndexStart:yIndexEnd,yIndexStart:yIndexEnd] = self.H_y[yIndexStart:yIndexEnd,yIndexStart:yIndexEnd] + amp

  def normalize_wavefunction(self):
    self.psi = (1/np.linalg.norm(self.psi))*self.psi
    
  def timeEvolution(self,tau,duration):
    self.duration = duration
    # Define A and B matrices
    A_x = sp.identity(self.sizeX) - tau/(1j)*self.H_x
    B_x = sp.identity(self.sizeX) + tau/(1j)*self.H_x
    # self.time_evolved_psi_x = np.zeros((self.sizeX,duration),dtype=complex)

    A_y = sp.identity(self.sizeY) - tau/(1j)*self.H_y
    B_y = sp.identity(self.sizeY) + tau/(1j)*self.H_y
    self.time_evolved_psi = np.zeros((self.sizeX,self.sizeY,duration),dtype=complex)

    # Start time evolution of particle
    for i in range(0,duration):
      print(i)
      # sys.stdout.write("Timestep: %d of %d \r" % (i+1,duration) )
      for ii in range(0,self.sizeY):
        # Solve linear equation Crank Nicelson
        # x-direction
        self.time_evolved_psi[:,ii,i],_ = linalg.bicgstab(A_x,B_x.dot(self.psi[:,ii]).transpose())
        self.psi[:,ii] = self.time_evolved_psi[:,ii,i]
        # y-direction
        self.time_evolved_psi[:,ii,i],_ = linalg.bicgstab(A_y,B_y.dot(self.psi[:,ii]).transpose())
        self.psi[:,ii] = self.time_evolved_psi[:,ii,i]

  def plot2D(self):
    time_evolved_probability = np.real(np.multiply(self.time_evolved_psi,np.conj(self.time_evolved_psi)))
    # time_evolved_probability_y = np.real(np.multiply(self.time_evolved_psi_y,np.conj(self.time_evolved_psi_y)))
    # time_evolved_probabilit y= np.zeros((self.sizeX,self.sizeY,self.duration),dtype=float)

    # for i in range(0,self.duration):
        # time_evolved_probability[:,:,i] = np.outer(time_evolved_probability_x[:,i],time_evolved_probability_y[:,i])

    # print(time_evolved_probability[:,:,0])

    # animate_wavefunction(time_evolved_probability,self.L,self.a,self.duration)
    # plt.plot(self.xAxis,sel)

    x_mesh,y_mesh = np.meshgrid(self.xAxis,self.yAxis)
    
    # fig = plt.figure()
    # ax = plt.axes(xlim=(0, self.L_x), ylim=(0, self.L_y))
    # plt.xlabel(r'x')
    # plt.ylabel(r'y')

    # anim = animation.FuncAnimation(fig, animate, frames=Nt)

    for i in range(0,self.duration):
        grid = time_evolved_probability[:,:,i]
        plt.imshow(grid, interpolation='none')
        plt.show()
