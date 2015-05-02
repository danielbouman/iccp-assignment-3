# ICCP assignment 3
import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

class Particle:
  k = 5.
  
  def __init__(self,a,L,sigma):
    self.L = L # box length
    self.a = a # lattice constant
    
    # Define the momentum operator matrix
    self.H = np.divide(sp.eye(a/L),-a**2)
    self.H += sp.eye(a/L,k=-1)
    self.H += sp.eye(a/L,k=1)
    
    # Wave function
    self.psi = 1/(sigma*np.sqrt(2*np.pi))*np.exp(i*self.k*x)
    
  def potential(pos,amp):
    # Add potential to hamiltonian matrix
    H_index = pos/self.a
    self.H(H_index,H_index) = amp
  
  def timeEvolution(tau,hbar):
    # Define A and B matrices
    self.A = sp.identity(a/L) - tau/(i*hbar)*self.H
    self.B = sp.identity(a/L) + tau/(i*hbar)*self.H
    
    # Time is run here
    for i in [1,2]:
      psi_new = linalg.bicgstab(self.A,self.B*psi)