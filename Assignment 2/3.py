from pylab import *
#Import the TDMA solver code
from tdma import tdma_solve
from time import sleep
#Define Grid Size
n = 20
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

#Line by line Gauss Seidel or Point by Point
lbl = 1

#Start the numerical solution of the problem
#Convergence criteria, Residual has to fall below Rtol
while Rmax > Rtol:
  
  #Intialise Residual for every iteration
  Rmax = 0
  #Check if we have to do line by line or not
  if lbl == 0:
    #Iterate over each point in Point by point Gauss Seidel
    for i in range(n-1,0,-1):
      for j in range(n-1,0,-1):
        #North, West, South, East Points
        Aw = k*dy/dx
        Ae = k*dy/dx
        An = k*dx/dy
        As = k*dx/dy
        #Current point index
        Ap = Aw + Ae + An + As
        #Generation term
        b = S*dx*dy
        #Calculate residual from previous time step and new values
        Rmax = max(abs(Aw*T[i][j-1] + Ae*T[i][j+1] + An*T[i+1][j] + 
               As*T[i-1][j] + b - Ap*T[i][j]) , Rmax)
        
        #Assign the new value
        T[i][j] = (Aw*T[i][j-1] + Ae*T[i][j+1] + An*T[i+1][j] + As*T[i-1][j] + b)/Ap
    
    #Boundary Conditions    
    for i in range(0,n+1):
      #Top and Right, temperature given
      T[i][n] = Ts
      T[n][i] = Ts
      #Left and Bottom, zero flux
      if (i < n) and (i > 0): 
        #Left face
        T[i][0] = (2*Ae*T[i][1] + An*T[i+1][0]/2 + As*T[i-1][0]/2 + S*dx*dy/2)/(2*Ae+An/2+As/2)
        #Bottom face
        T[0][i] = (2*An*T[1][i] + Aw*T[0][i-1]/2 + Ae*T[0][i+1]/2 + S*dx*dy/2)/(2*An+Aw/2+Ae/2)
    #Left bottom point
    T[0][0] = (An*T[1][0] + Ae*T[0][1] + S*dx*dy/4)/(An+Ae)
    
  else:   
    #Line By Line Gauss Seidel
    #Top and Line Boundary conditions
    for i in range(0,n+1):
      T[n][i] = Ts  
      
    #Iterate over every line  
    for i in range(n-1,-1,-1):
      #Assign indices for internal grid points
      for j in range(1,n):
        #West and East points for TDMA
        Ad[j] = -k*dy/dx
        Au[j] = -k*dy/dx
        #North and South points come in the constant
        An = k*dx/dy
        As = k*dx/dy
        #Bottom line TDMA
        if i==0:
          Ad[j] /= 2
          Au[j] /= 2
          An *= 2
          As = 0
        
        #Current point index
        A[j] = An+As-Ad[j]-Au[j]
        #Constant for TDMA
        if i==0:
          b[j] = S*dx*dy/2 + An*T[i+1][j]
        else:
          b[j] = S*dx*dy + An*T[i+1][j] + As*T[i-1][j]
        #Calculate Residual from previous time step values
        if i!=0:
          Rmax = max(abs(A[j]*T[i][j] + Ad[j]*T[i][j-1] + Au[j]*T[i][j+1]-b[j]), Rmax)
      
      #Right Boundary point
      b[n] = Ts
      A[n] = 1
      #Left Boundary point
      if i==0:
        Au[0] = -k*dy/dx
        A[0] = -Au[0]+An
        b[0] = S*dx*dy/4 + An*T[1][0]
      else:
        Au[0] = -2*k*dy/dx
        A[0] = -Au[0]+An/2+As/2
        b[0] = S*dx*dy/2 + An*T[i+1][0]/2 + As*T[i-1][0]/2
      
      #Solve using TDMA
      TD = tdma_solve((A, Ad, Au), b)
      #Copy the result into the temperature matrix
      for a in range(0,n+1):
        T[i][a] = TD[a]
  
  #Solve the matrix and keep track on the residuals
  iter += 1
  print iter, Rmax

#Analytical solution, for all points
for i in range(0,n+1):
  for j in range(0,n+1):
    #Series constants
    sum = 0
    delta = 1
    N = 1
    #Analytical solution's series convergence criteria
    while abs(delta) > Rtol:
      la = (2*N-1)*pi/2
      #Series term
      delta = cos(la*x[j])*cosh(la*y[i])/((la**3)*cosh(la))
      sum += ((-1)**N)*delta
      N += 1
    #Calculate T
    #T[i][j] = Ts + S*(1-(x[j])**2)/(2*k) + sum*2*S/k
  
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


  

