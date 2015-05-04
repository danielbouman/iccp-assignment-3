import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

class Particle:
  
  def __init__(self,a,L,delta,k):
    self.L = L # box length
    self.a = a # spacial resolution
    
    # Define the momentum operator matrix
    self.H = np.divide(sp.eye(L/a),-a**2)
    self.H += sp.eye(L/a,k=-1)
    self.H += sp.eye(L/a,k=1)
    
    # Wave function
    self.psi = 1/(sigma*np.sqrt(2*np.pi))*np.exp(i*self.k*x)
    
  def potential(pos,amp):
    # Add potential to hamiltonian matrix
    H_index = pos/self.a
    self.H(H_index,H_index) = amp
  
  def timeEvolution(tau,hbar,duration):
    # Define A and B matrices
    A = sp.identity(self.L/self.a) - tau/(i*hbar)*self.H
    B = sp.identity(self.L/self.a) + tau/(i*hbar)*self.H
    
    # Time is run here
    for i in range(0,duration+1):
      psi_new = linalg.bicgstab(A,B*psi)
    