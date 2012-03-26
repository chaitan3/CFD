from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
from time import sleep
#Define Grid Size
n = 10
#Grid spacing, unit: m
l = 0.1
dx = l/n
#Conductivity of the material, unit: W/m-K
kc = 0.05
beta = 0.001
#Heat Transfer coefficients on the left and right faces, unit: W/m2-K
h_r = 5.0
h_l = 30.0
#Temperature of the surroundings of left and right faces, unit: K
T_l = 300.0
T_r = 30.0
#Vector of the x-coordinate nodes
x = arange(0, 0.1 + dx, dx)

#Function to calculate the conductivity of the material given the temperature
def k(T):
  return kc*(1+beta*T)

#Solution vector initialisation
T = 300*ones((n+1,1))
A = zeros((n+1,1))
Au = ones((n+1,1))
Ad = ones((n+1,1))

#Residual tolerance
Rtol = 1e-6
Rmax = 1
cont = True

while Rmax > Rtol:
  Rmax = 0
  #Calculate the Conductivities based on the T calculated in the previous step
  ks = k(T)
  #Filling up the diagonal of the matrix to be solved
  for i in range(1,n):
    #Interface conductivities: harmonic mean
    kw = 2*ks[i]*ks[i-1]/(ks[i]+ks[i-1])
    ke = 2*ks[i]*ks[i+1]/(ks[i]+ks[i+1])
    A[i] = (kw+ke)/dx
    Ad[i] = -kw/dx
    Au[i] = -ke/dx
    #Calculate point residuals and get the maximum
    R = abs(A[i]*T[i] + Ad[i]*T[i-1] + Au[i]*T[i+1])
    Rmax = max(R, Rmax)
  
  #Boundary Conditions: Heat Flux
  A[0] = ks[0]/dx + h_l
  A[n] = ks[n]/dx + h_r
  #Left Boundary condition
  Au[0] = -ks[0]/dx
  #Right Boundary condition
  Ad[n] = -ks[n]/dx
  
  #Source term
  b = zeros((n+1,1))
  #Boundary conditions
  b[0] = h_l*T_l
  b[n] = h_r*T_r
  
  #Residuals for boundary points
  R = abs(A[n]*T[n] + Ad[n]*T[n-1] - b[n])
  Rmax = max(R, Rmax)
  R = abs(A[0]*T[0] + Au[0]*T[1] - b[0])
  Rmax = max(R, Rmax)
  
  #Solve the matrix and keep track on the residuals
  print Rmax
  T = tdma_solve((A, Ad, Au), b)
  
#Plot the numerical solution
plot(x,T)
xlabel('x(m)')
ylabel('T(degrees Celcius)')
show()


  
