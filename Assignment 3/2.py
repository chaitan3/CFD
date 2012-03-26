from pylab import *
from time import sleep
#Define Grid Size
n = 20
#Grid spacing, dx, dy: 2D system
l = n*1.0
dx = 1.0
dy = 1.0
#Time step
dt = 0.5
#Density and diffusivity of the material
rho = 1.0
T = 1.0
#Source term constants
a = 10.0
b = 2.0
#Convection in x and y direction
u = 1.0
v = 4.0 
Fx = rho*u
Fy = rho*v
#phi constants at the walls
phi_t = 0.0
phi_0 = 100.0
#Vector of the x,y-coordinate nodes
x = arange(0, l + dx, dx)
y = arange(0, l + dy, dy)

#Variable to choose between the schemes
#f = 0, CDS
#f = 1, UDS
#f = 2, Hybrid Scheme
f = 2

#Intial time step value for phi
phi = 50.0*ones((n+1,n+1))

#Residual tolerance and Convergence criteria for steady state
Rtol = 1e-3
#Convergence criteria variable for steady state
tRmax = 1

#Inialise time
t = 0

#Run Steady state or transient
steady = 0

#Generate the mesh for the plot
(X,Y) = meshgrid(x,y)

#Boundary Conditions    
for i in range(0,n+1):
  #phi specified at the walls, no need for this to be inside the iteration loop
  phi[i][n] = phi_t
  phi[n][i] = phi_t
  phi[0][i] = phi_0
  phi[i][0] = phi_0

#Start the numerical solution of the problem
#Overall Convergence criteria for steady state
while tRmax > Rtol:
  #Print the current time step
  print t, tRmax
  #Copy the old time step phi
  phio = phi.copy()
  #Plot the numerical solution
  
  #Contour plot
  #contourf(X,Y, phi, arange(phi_t,phi_0+Rtol,10))
  #xlabel('x(m)')
  #ylabel('y(m)')
  #colorbar()
  #savefig('2t'+str(t)+'.pdf')
  #clf()
  
  #Initialise the convergence criteria and residuals
  tRmax = 0
  Rmax = 1
  #Iteration number
  iter = 0
  
  #Residual has to fall below Rtol for convergence of one time step
  while Rmax > Rtol:
    #Intialise Residual for every iteration
    Rmax = 0
    #Iterate over each point in Point by point Gauss Seidel
    for i in range(1,n):
      for j in range(1,n):
        #North, West, South, East Points
        if f == 0:
          #CDS
          Aw = (T/dx+Fx/2)*dy
          Ae = (T/dx-Fx/2)*dy
          An = (T/dy-Fy/2)*dx
          As = (T/dy+Fy/2)*dx
        else:
          #UDS
          if f == 1:
            Aw = (T/dx+Fx)*dy
            Ae = (T/dx)*dy
            An = (T/dy)*dx
            As = (T/dy+Fy)*dx
          else:
            #Hybrid
            Aw = max(T/dx+Fx/2,Fx)*dy
            Ae = max(T/dx-Fx/2,0)*dy
            An = max(T/dy-Fy/2,0)*dx
            As = max(T/dy+Fy/2,Fy)*dx
          
        #Current point index
        Ap = Aw + Ae + An + As + b*dx*dy 
        #Generation term
        bc = a*dx*dy 
        if steady == 0:
          Ap += rho*dx*dy/dt
          #Generation term it has the previous time step value(implicit scheme)
          bc += rho*phio[i][j]*dx*dy/dt
          
        #Calculate current residual using values from current and previous iteration depending on the sweep
        tmp = Aw*phi[i][j-1] + Ae*phi[i][j+1] + An*phi[i+1][j] + As*phi[i-1][j] + bc - Ap*phi[i][j]
        Rmax = max(abs(tmp) , Rmax)
        
        #Assign the new value
        phi[i][j] += 1*(tmp/Ap)
        #update steady state convergence
        tRmax = max(abs(phi[i][j]-phio[i][j])/dt, tRmax)
    
    #Solve the matrix and keep track on the residuals
    iter += 1
    print iter, Rmax
    
  #Increment time by the time step
  t += dt
  
  #Check if Steady State simulation
  if steady:
    break

#Plot the numerical solution
#Contour plot
contourf(X,Y, phi, arange(phi_t,phi_0+Rtol,10))
colorbar()
xlabel('x(m)')
ylabel('y(m)')
#savefig('2t'+str(t)+'s.pdf')
show()


  

