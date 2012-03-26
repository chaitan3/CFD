from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
from time import sleep
from scipy import weave
#Define Grid Size
n = 100
#Grid spacing, dx, dy: 2D system
l = 1.0
dx = l/n
dy = l/n
#Conductivity of the material
k = 0.5
#Source Term, unit: W/m2
S = 1000
#Temperature of the surroundings of left and right faces, unit: K
Ts = 100
#Vector of the x,y-coordinate nodes
x = arange(0, l + dx, dx)
y = arange(0, l + dy, dy)

#Intial Guess value for temperature
T = 5*Ts*ones((n+1,n+1))
#Vectors for TDMA
A = zeros((n+1,1))
Au = zeros((n+1,1))
Ad = zeros((n+1,1))
b = zeros((n+1,1))

#Residual tolerance
Rtol = 1e-3
#Residual maxima
Rmax = 1

#Iteration number
iter = 0

#Start the numerical solution of the problem
#Convergence criteria, Residual has to fall below Rtol
R = T[1:-1,1:-1].copy()

while Rmax > Rtol:
  
  #Intialise Residual for every iteration
  Aw = k*dy/dx
  Ae = k*dy/dx
  An = k*dx/dy
  As = k*dx/dy
  Ap = Aw + Ae + An + As
  b = S*dx*dy
    
  e = "T[1:-1,1:-1] = (Aw*T[1:-1,:-2] + Ae*T[1:-1,2:] + An*T[2:,1:-1] + As*T[:-2,1:-1] + b)/Ap"
  weave.blitz(e, check_size=0)
  e = "R = Aw*T[1:-1,:-2] + Ae*T[1:-1,2:] + An*T[2:,1:-1] + As*T[:-2,1:-1] + b - Ap*T[1:-1,1:-1]"
  weave.blitz(e, check_size=0)
  e = "T[:,n] = Ts"
  weave.blitz(e, check_size=0)
  e = "T[n,:] = Ts"
  weave.blitz(e, check_size=0)
  e = "T[:,0] = T[:,1]"
  weave.blitz(e, check_size=0)
  e = "T[0,:] = T[1,:]"
  weave.blitz(e, check_size=0)
  
  R = abs(R)
  Rmax = R.max()

  #Solve the matrix and keep track on the residuals
  iter += 1
  if iter % 100 == 0:
    print iter, Rmax
print iter, Rmax

#Plot the numerical solution
#Generate the mesh
(X,Y) = meshgrid(x,y)
#Upper limit
Tsmax = 700
T[n][n]= Ts + Rtol
#Contour plot
contourf(X,Y, T,arange(Ts,Tsmax,5))
#Set the colour bar
colorbar(ticks=arange(Ts,Tsmax,100))
#Axis limits
axis([0,1,0,1])
xlabel('x(m)')
ylabel('y(m)')
#show()


  

