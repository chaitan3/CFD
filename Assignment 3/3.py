from pylab import *
from time import sleep
#Grid spacing, dx, dy: 2D system
l = 20.0
dx = 0.1
dy = 0.1
#Height of channel
h = 1.0
#Number of grid points in each direction
nx = int(l/dx)
ny = int(h/dy)
#Material Properties: Density, specific heat, conductivity 
rho = 1.0
c = 100.0
k = 1.0
#Wall Flux and inlet temperature
hf = 1000.0
Ti = 50.0
#Vector of the x,y-coordinate nodes
x = arange(0, l + dx, dx)
y = arange(0, h + dy, dy)

#Intial Guess value for temperature
T = Ti*ones((ny+1,nx+1))

#Horizontal velocity function
def velocity(yp):
  yp -= h/2
  return 1.5*(1-4*yp*yp)

#Residual tolerance
Rtol = 1e-3
#Residual maxima
Rmax = 1

#Function to update the current grid point value in Gauss Seidel
def update(Aw, Ae, An, As, b, i, j):
  global Rmax, T, nx, ny
  #Grid point coefficient
  Ap = Aw + Ae + An + As
  #constant term
  tmp = b
  #add all the neighbouring point terms provided they exist
  if  abs(Aw) > Rtol:
    tmp += Aw*T[i][j-1]
  if  abs(Ae) > Rtol:
    tmp += Ae*T[i][j+1]
  if  abs(An) > Rtol:
    tmp += An*T[i+1][j]
  if  abs(As) > Rtol:
    tmp += As*T[i-1][j]
    
  tmp -= Ap*T[i][j]
  #Find max residual
  Rmax = max(abs(tmp) , Rmax)
  #Assign The new value
  T[i][j] += 1*(tmp/Ap)

#Iteration number
iter = 0

#Start the numerical solution of the problem
#Convergence criteria, Residual has to fall below Rtol
while Rmax > Rtol:
  #Intialise Residual for every iteration
  Rmax = 0
  #Iterate over each point in Point by point Gauss Seidel
  for i in range(1,ny):
    for j in range(1,nx):
      #North, West, South, East Points
      #Hybrid scheme used
      #Assign flow rate in horizontal direction
      F = rho*c*velocity(y[i])
      #Only x-direction has flow, hence only W and E have flow in their coefficients
      Aw = max(F, k/dx+F/2, 0)*dy
      Ae = max(-F, k/dx-F/2, 0)*dy
      An = k*dx/dy
      As = k*dx/dy
      #Update the value at the current grid point
      update(Aw, Ae, An, As, 0, i, j)
  
  #Boundary conditions
  #bottom left point, inlet
  T[0][0] = Ti
  #top left point, inlet
  T[ny][0] = Ti
  #bottom right point, fully developed flow BC(adiabatic dT/dx = 0)
  Aw = k*dy/(2*dx)
  An = k*dx/(2*dy)
  As = 0
  Ae = 0
  b = hf*dx/2
  update(Aw, Ae, An, As, b, 0, nx)
  
  #top right point, fully developed flow
  Aw = k*dy/(2*dx)
  An = 0
  As = k*dx/(2*dy)
  Ae = 0
  b = hf*dx/2
  update(Aw, Ae, An, As, b, ny, nx)
  
  for i in range(1,ny):
    #inlet temperature
    T[i][0] = Ti
    
    #outlet, fully developed flow
    F = rho*c*velocity(y[i])
    Aw = max(F, k/dx+F/2, 0)*dy
    An = k*dx/(2*dy)
    As = k*dx/(2*dy)
    Ae = 0
    b = 0
    update(Aw, Ae, An, As, b, i, nx)
  
  for i in range(1,nx):
    #bottom wall, heat flux given
    Aw = k*dy/(2*dx)
    Ae = k*dy/(2*dx)
    An = k*dx/dy
    As = 0
    b = hf*dx
    update(Aw, Ae, An, As, b, 0, i)
    
    #top wall, heat flux given
    Aw = k*dy/(2*dx)
    Ae = k*dy/(2*dx)
    As = k*dx/dy
    An = 0
    b = hf*dx
    update(Aw, Ae, An, As, b, ny, i)
        
  #Solve the matrix and keep track on the residuals
  iter += 1
  print iter, Rmax

#Plot the numerical solution
#Generate the mesh
(X,Y) = meshgrid(x,y)
#Contour plot
contourf(X,Y, T, arange(Ti,750,1))
colorbar()
xlabel('x(m)')
ylabel('y(m)')
axis([0, l, 0, h])

#Initialise Bulk temperature, heat transfer coefficient and nusselt number arrays
Tb = x.copy()
ht = x.copy()
Nu = x.copy()
for i in range(0,nx+1):
  #Calculate the bulk temperature
  Tb[i] = 0
  #Simpsons Integration over all y for each x
  for j in range(1,ny):
    if j%2 == 0:
      f = 2
    else:
      f = 4
    Tb[i] += f*T[j][i]*velocity(y[j])
  #Bulk temperature
  Tb[i] *= dy/3
    
  #Heat flux calculation for heat transfer coefficient
  #Wall heat flux, using the bottom wall
  Jn = -k*(T[1][i] - T[0][i])/dy
  Je = 0
  Jw = 0
  #left and right points
  if i != nx:
    Je = -k*(T[0][i+1] - T[0][i])/dx
  if i != 0:
    Jw = -k*(T[0][i] - T[0][i-1])/dx
  #Wall heat flux
  Js = Jn + (Je-Jw)*dy/(2*dx)
  if i == nx:
    Js = Jn - Jw*dy/dx
  if i == 0:
    Js = 0
  #Heat transfer coefficient
  ht[i] = Js/(T[0][i] - Tb[i])
  #Nusselt number
  Nu[i] = ht[i]*2*h/k
#plot(x, Tb)
#plot(x, ht)
#plot(x, Nu)
#plot(x, 8.23*ones(nx+1))
#print Nu
#show()


  


