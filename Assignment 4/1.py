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
#Inlet velocity
u0 = 1.0
mu = 0.01

ua = 0.5
pa = 0.8


#Vector of the x,y-coordinate nodes
x = arange(0, l + dx, dx)
y = arange(0, h + dy, dy)

#Intial Guess value for temperature
u = u0*ones((ny+1,nx+1))
bm = zeros((ny+1,nx+1))
p = zeros((ny+1,nx+1))
d = zeros((ny+1,nx+1))

for j in range(0, nx+1):
  u[0][j] = 0
  u[ny][j] = 0

#Mass source tolerance
btol = 1e-5
#Mass source
bmax = 1


#Start the numerical solution of the problem
#Convergence criteria, mass source has to fall below btol
while 1:
  for iter in range(0,4):
    bmax = 0
    #Iterate over each point in Point by point Gauss Seidel
    
    for i in range(1,ny):
      for j in range(1,nx):
        #North, West, South, East Points
        #Hybrid scheme used
        #Assign flow rate in horizontal direction
        
        #Only x-direction has flow, hence only W and E have flow in their coefficients
        if j == 0:
          Aw = 2*mu*dy/dx
        else:
          F = rho*(u[i][j-1]+u[i][j])/2
          Aw = max(F, mu/dx+F/2, 0)*dy
        if j == nx-1:
          Ae = 0
        else:
          F = rho*(u[i][j]+u[i][j+1])/2
          Ae = max(-F, mu/dx-F/2, 0)*dy
        
        An = mu*dx/dy
        As = mu*dx/dy
        Ap = Aw + Ae + An + As
        b = (p[i][j]-p[i][j+1])*dy
        if j == 0:
          b += (rho*u0*u0 + 2*mu*u0/dx)*dy
          Aw = 0
        #Update the value at the current grid point
        
        u[i][j] = (b + Aw*u[i][j-1] + Ae*u[i][j+1] + An*u[i+1][j] + As*u[i-1][j] + (1-ua)*Ap*u[i][j]/ua)*ua/Ap
        bm[i][j] = rho*(u[i][j-1]-u[i][j])*dy
        bmax = max(bm[i][j], bmax)
        d[i][j] = dy/Ap
          
    #Solve the matrix and keep track on the residuals
    print 'u',iter, bmax
    
  if bmax < btol:
    break
    
  pp = zeros((ny+1,nx+1))
  for iter in range(0, 10):
    for i in range(1,ny):
      for j in range(1,nx):
        Ae = rho*d[i][j]*dy
        Aw = rho*d[i][j-1]*dy
        b = bm[i][j]
        #Update the value at the current grid point
        Ap = Aw + Ae
        pp[i][j] = (Aw*pp[i][j-1]+Ae*pp[i][j+1]+b)/Ap
    for j in range(1, nx):
      pp[0][j] = pp[1][j]
      pp[ny][j] = pp[1][j]
    for i in range(0,ny+1):
      pp[i][0] = pp[i][1]
      
      Aw = rho*d[i][nx-1]*dy
      b = bm[i][nx]
      #Update the value at the current grid point
      Ap = Aw
      
      pp[i][nx] = (Aw*pp[i][nx-1]+b)/Ap
      
    print 'pp', iter    
  for i in range(1,ny):
    for j in range(1,nx):
      u[i][j] += d[i][j]*(pp[i][j]-pp[i][j+1])
      p[i][j] += pa*pp[i][j]
  
  
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


  



