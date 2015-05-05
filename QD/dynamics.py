import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from anim2D import animate_wavefunction
import pylab
from mpl_toolkits.mplot3d import Axes3D

class Particle:
  
  def __init__(self,a,L_x,L_y,sigma_x,sigma_y,k_x,k_y,mu_x,mu_y):
    self.L_x = L_x # box length
    self.L_y = L_y # box length
    self.a = a # spatial resolution
    self.xAxis = np.linspace(0,L_x,L_x/a)
    self.yAxis = np.linspace(0,L_y,L_y/a)
    
    # Define the momentum operator matrices
    self.H_x = np.divide(sp.diags([1,-2,1],[-1,0,1], shape=(L_x/a,L_x/a)),-a**2).todense()
    self.H_y = np.divide(sp.diags([1,-2,1],[-1,0,1], shape=(L_y/a,L_y/a)),-a**2).todense()

    # Wave function
    self.psi_x = np.multiply(np.exp((-1/sigma_x)*np.power(self.xAxis-mu_x,2)),np.exp(-1j*k_x*self.xAxis))
    self.psi_y = np.multiply(np.exp((-1/sigma_y)*np.power(self.yAxis-mu_y,2)),np.exp(-1j*k_y*self.yAxis))
    
  def potential(self,xStart,xEnd,yStart,yEnd,amp):
    # Add potential to hamiltonian matrix
    xIndexStart = xStart/self.a
    xIndexEnd = xEnd/self.a
    yIndexStart = yStart/self.a
    yIndexEnd = yEnd/self.a
    self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] = self.H_x[xIndexStart:xIndexEnd,xIndexStart:xIndexEnd] + amp
    self.H_y[yIndexStart:yIndexEnd,yIndexStart:yIndexEnd] = self.H_y[yIndexStart:yIndexEnd,yIndexStart:yIndexEnd] + amp

  def normalize_wavefunction(self):
    inner_product = self.a*sum(np.multiply(self.psi_x,np.conj(self.psi_x)))*self.a*sum(np.multiply(self.psi_y,np.conj(self.psi_y)))
    self.psi_x = (1/np.sqrt(np.abs(inner_product)))*self.psi_x
    
  def timeEvolution(self,tau,hbar,duration):
    self.duration = duration
    # Define A and B matrices
    A_x = sp.identity(self.L_x/self.a) - tau/(1j*hbar)*self.H_x
    B_x = sp.identity(self.L_x/self.a) + tau/(1j*hbar)*self.H_x
    self.time_evolved_psi_x = np.zeros((self.L_x/self.a,duration),dtype=complex)

    A_y = sp.identity(self.L_y/self.a) - tau/(1j*hbar)*self.H_y
    B_y = sp.identity(self.L_y/self.a) + tau/(1j*hbar)*self.H_y
    self.time_evolved_psi_y = np.zeros((self.L_x/self.a,duration),dtype=complex)

    # Time is run here
    for i in range(0,duration):
        # Solve linear equations for each dimension
        self.time_evolved_psi_x[:,i],_ = linalg.bicgstab(A_x,B_x.dot(self.psi_x).transpose())
        self.psi_x = self.time_evolved_psi_x[:,i]
        self.time_evolved_psi_y[:,i],_ = linalg.bicgstab(A_y,B_y.dot(self.psi_y).transpose())
        self.psi_y = self.time_evolved_psi_y[:,i]

  def plot2D(self):
    time_evolved_probability_x = np.real(np.multiply(self.time_evolved_psi_x,np.conj(self.time_evolved_psi_x)))
    time_evolved_probability_y = np.real(np.multiply(self.time_evolved_psi_y,np.conj(self.time_evolved_psi_y)))
    time_evolved_probability = np.zeros((self.L_x/self.a,self.L_y/self.a,self.duration),dtype=float)

    for i in range(0,self.duration):
        time_evolved_probability[:,:,i] = np.outer(time_evolved_probability_x[:,i],time_evolved_probability_y[:,i])

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
