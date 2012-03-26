from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
#Define Grid Size
n = 10
#Grid spacing
dx = 0.1/n
#Conductivity of the material
k = 0.1
#Heat Transfer coefficients on the left and right faces
h_r = 5.0
h_l = 30.0
#Temperature of the surroundings of left and right faces
T_l = 300.0
T_r = 30.0
#Vector of the x-coordinate nodes
x = arange(0, 0.1 + dx, dx)

#Filling up the diagonal of the matrix to be solved
A = 2*ones((n+1, 1))
#Boundary Conditions: Heat Flux
A[0] = k/dx + h_l
A[n] = k/dx + h_r
#Filling up the upper diagonal of the matrix to be solved,
#constant conductivity
Au = -1*ones((n+1, 1))
#Left Boundary condition
Au[0] = -k/dx
#Filling up the lower diagonal of the matrix to be solved
Ad = -1*ones((n+1, 1))
#Right Boundary condition
Ad[n] = -k/dx
#The constant terms of the Linear equation to be solved, 
#no source term
b = zeros((n+1,1))
#Boundary conditions
b[0] = h_l*T_l
b[n] = h_r*T_r
#Solution vector
T = zeros((n+1,1))

#Code for plotting the analytical solution
c2 = 292.70270270270271
c1 = -300*(300-c2)
xt = arange(0, 0.1 + 0.1/200, 0.1/200)
y=c1*xt+c2
#Plot the analytical solution
#plot(xt, y, color='red')
#Solve the linear equation
T = tdma_solve((A, Ad, Au), b)
#Plot the numerical solution
plot(x,T)
xlabel('x(m)')
ylabel('T(degrees Celcius)')
#The Boundary values of numerical and analytical solutions
print T[0], T[n]
print y[0], y[-1]
show()
