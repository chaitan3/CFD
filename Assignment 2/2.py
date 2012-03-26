from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
from time import sleep
#Define Grid Size
n = 10
l = 1.0
#Grid spacing, Time step, unit: m, s
dx = l/n
dt = dx*dx/2
#dt = 0.005
#Thermal Conductivity of the material, unit: m2/hr
alpha = 1.0
#Temperature of the surroundings of left and right faces, unit: K
T_r = 300.0
T_l = 100.0
#Choosing Explicit/Crank-Nicholson/Implicit Schemes by changing factor f
f = 1
#Vector of the x-coordinate nodes
x = arange(0, l + dx, dx)

print dx, dt

#Initial T values at time = 0
T = 100+200*(x**2)
#Initialise the time variable to 0
t = 0
#Relative tolerance for the solution
Rtol = 1e-3
#Tolerance check variable
Rmax = 1

#TDMA vector initialisation
A = zeros((n+1,1))
Au = zeros((n+1,1))
Ad = zeros((n+1,1))
b = zeros((n+1,1))

#Function for calculating Analytical solution
def analytic(t):
  #Get the global values of the various problem constants
  global T_r, T_l, alpha, x
  #Declare the temperature array
  T = x.copy()
  #Calculate T for every x point
  for i in range(0,len(x)):
    #Constants for the sum
    su = 0
    delta = 1
    N = 1
    #Analytical solution's series convergence criteria
    while abs(delta) > Rtol:
      la = (2*N-1)*pi/2
      #Series term
      delta = cos(la*x[i])*(e**(-(la**2)*alpha*t))/(la**3)
      su += ((-1)**(N+1))*delta
      N += 1
    #Calculate T
    T[i] = T_r + 4*(T_l-T_r)*su
  #Return the calculated solution
  return T

#Start the numerical solution of the problem
#Convergence criteria, Residual has to fall below Rtol
while Rmax > Rtol:
  
  #Copy the previous time step solution for future use
  To = T.copy()
  #If f = 0, explicit solution(no TDMA required), f > 0, as new values are used TDMA required.
  if f == 0:
    #Iterate over all inside grid points
    for i in range(1,n):
      #Calculate T from previous time step values
      T[i] = (dt/dx)*(To[i-1]*alpha/dx + To[i+1]*alpha/dx + (dx/dt - 2*alpha/dx)*To[i])
    #Boundary point no flux Half CV method
    T[0] = (dt/dx)*(To[1]*alpha/dx + (dx/dt - alpha/dx)*To[0])
    #Calculate the difference between successive time step solutions by taking norm
    Rmax = norm(T-To)
  else:
    #Filling up the diagonal of the matrix to be solved
    for i in range(1,n):
      A[i] = 2*f*alpha/dx+dx/dt
      Ad[i] = -f*alpha/dx
      Au[i] = -f*alpha/dx
      b[i] = dx*T[i]/dt+(1-f)*alpha*(T[i-1]+T[i+1]-2*T[i])/dx
    #Right Boundary Condition: Temperature specified
    A[n] = 1
    b[n] = T_r
    #Left Boundary Conditions: Heat Flux zero
    Au[0] = -f*alpha/dx
    A[0] = f*alpha/dx+dx/dt
    b[0] = dx*T[0]/dt+(1-f)*alpha*(T[1]-T[0])/dx
    #Au[0] = -1
    #A[0] = 1
    
    #Solve using TDMA
    T = tdma_solve((A, Ad, Au), b)
    #Calculate the difference between successive time step solutions by taking norm
    Rmax = norm(T-To)
  
  #Increment t by the time step
  t += dt
  print t,Rmax
  
  #Plot the solution every 0.25 second
  if abs(4*t-round(t*4)) < Rtol:
    #Plot numerical solution
    plot(x,T, 'r--')
    #Calculate analytic solution
    Ta = analytic(t)
    #Check the two solutions
    print T[0], Ta[0]
    #Plot Analytical solution
    plot(x,Ta)
    #show()
  
#Plot the numerical solution
plot(x,T)
xlabel('x(m)')
ylabel('T(degrees Celcius)')
axis([0,1,180,T_r])
show()


  

