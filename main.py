# ICCP assignment 3
import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as linalg

class Particle:
  k = 5.
  
  def __init__(self,a,L):
    self.L = L # box length
    self.a = a # lattice constant
    
    # Define the momentum operator matrix
    self.H = sp.eye(a/L)
    self.H += sp.eye(a/L,k=-1)
    self.H += sp.eye(a/L,k=1)
    
  def wavefunction(sigma,mu):
    self.psi = 1/(sigma*np.sqrt(2*np.pi))*np.exp(i*k*x)
    
  def potential(x):
    