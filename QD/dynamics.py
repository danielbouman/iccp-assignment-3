import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

class Particle:
  
  def __init__(self,a,L,sigma,k):
    self.L = L # box length
    self.a = a # spatial resolution
    xAxis = np.linspace(0,L,L/a+1)
    
    # Define the momentum operator matrix
    self.H = np.multiply(sp.eye(L/a,k=0),-2)
    self.H = self.H + sp.eye(L/a,k=-1)
    self.H = self.H + sp.eye(L/a,k=1)
    self.H = np.divide(self.H,-a**2).todense()
    
    # Wave function
    self.psi = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-np.power(xAxis,2))
    
  def potential(pos,amp):
    # Add potential to hamiltonian matrix
    H_index = pos/self.a
    self.H[H_index,H_index] = amp
  
  def timeEvolution(tau,hbar,duration):
    # Define A and B matrices
    A = sp.identity(self.L/self.a) - tau/(i*hbar)*self.H
    B = sp.identity(self.L/self.a) + tau/(i*hbar)*self.H
    
    time_evolved_psi = np.zeros(L/a,duration,dtype=float)

    # Time is run here
    for ii in range(0,duration):
      psi_new = linalg.bicgstab(A,B*psi)
      time_evolved_psi[:,ii] = psi_new
    
