from pylab import *
#Function implementing the TDMA algorithm
def tdma_solve((A, Ad, Au), b):
  #Get size of linear equations to be solved
  n = size(b) - 1
  #Initialise the solution
  T = zeros((n+1,1))
  #Forward Elimination
  for i in range(1, n+1):
    f = Ad[i]/A[i-1]
    A[i] -= f*Au[i-1]
    b[i] -= f*b[i-1]
  #Back substitution
  T[n] = b[n]/A[n]
  for i in range(n-1,-1, -1):
    T[i] = (b[i]-Au[i]*T[i+1])/A[i]
  return T
