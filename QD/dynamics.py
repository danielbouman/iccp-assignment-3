import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from custom_functions import frange

class Particle:
  
  def __init__(self,a,L,sigma,k):
    self.L = L # box length
    self.a = a # spatial resolution
    xAxis = frange(0,L,a)
    
    # Define the momentum operator matrix
    self.H = np.multiply(sp.eye(L/a),-2)
    self.H += sp.eye(L/a,k=-1)
    self.H += sp.eye(L/a,k=1)
    self.H = np.divide(self.H,-a**2).todense()
    
    # Wave function
    self.psi = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-k*xAxis**2)
    
  def potential(pos,amp):
    # Add potential to hamiltonian matrix
    H_index = pos/self.a
    self.H[H_index,H_index] = amp
  
  def timeEvolution(tau,hbar,duration):
    # Define A and B matrices
    A = sp.identity(self.L/self.a) - tau/(i*hbar)*self.H
    B = sp.identity(self.L/self.a) + tau/(i*hbar)*self.H
    
    # Time is run here
    for i in range(0,duration+1):
      psi_new = linalg.bicgstab(A,B*psi)
    
