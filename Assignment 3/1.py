from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
#Define Grid Size
n = 400
#Grid spacing
l = 1.0
dx = l/n
#Density and diffusivity of the material
rho = 1.0
T = 0.05
#Convection: velocity and flow rate
u = 2
F = rho*u
#Phi values of left and right faces
phi_0 = 0.0
phi_l = 1.0
#Vector of the x-coordinate nodes
x = arange(0, l + dx, dx)

#Variable to decide whether to use CDS or UDS
uds = 1

#UDS
if uds == 1:
  #Filling up the diagonal of the matrix to be solved
  A = (2*(T/dx)+F)*ones((n+1, 1))
  #Filling up the upper diagonal of the matrix to be solved,
  Au = -1*(T/dx)*ones((n+1, 1))
  #Filling up the lower diagonal of the matrix to be solved
  #In upwind difference scheme only one of the surrounding values is used for convection
  Ad = -1*(T/dx+F)*ones((n+1, 1))
#CDS
else:
  #Filling up the diagonal of the matrix to be solved
  A = 2*(T/dx)*ones((n+1, 1))
  #Filling up the upper diagonal of the matrix to be solved,
  #In central difference scheme both of the surrounding values are used for convection
  Au = -1*(T/dx-F/2)*ones((n+1, 1))
  #Filling up the lower diagonal of the matrix to be solved
  Ad = -1*(T/dx+F/2)*ones((n+1, 1))
#The constant terms of the Linear equation to be solved, 
#no source term
b = zeros((n+1,1))
#Boundary conditions
#Left Boundary condition
Au[0] = 0
b[0] = phi_0
A[0] = 1
#Right Boundary condition
Ad[n] = 0
b[n] = phi_l
A[n] = 1
#Solution vector
phi = zeros((n+1,1))

#Code for plotting the analytical solution
P = rho*u*l/T
y = []
xt = arange(0, l + l/200, 0.1/200)
for i in xt:
  y.append(phi_0+(phi_l-phi_0)*(e**(P*i/l)-1)/(e**P-1))
#Plot the analytical solution
plot(xt, y, color='red')
#Solve the linear equation
phi = tdma_solve((A, Ad, Au), b)
#Plot the numerical solution
plot(x,phi)
xlabel('x(m)')
ylabel('phi')
axis([0, 1, 0, 1])
#show()

