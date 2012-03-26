from pylab import *
from tdma import tdma_solve
n = 100
l = 0.02
dx = l/n
k = 0.5
#Source term
S = 1e6
x = arange(0, l + dx, dx)

A =2*(k/dx)*ones((n+1, 1))
#Boundary conditions: Temperature given
A[0] = 1
A[n] = 1
Au = -(k/dx)*ones((n+1, 1))
Au[0] = 0
Ad = -(k/dx)*ones((n+1, 1))
Ad[n] = 0
#Constant term of linear equation, source term
b = (S*dx)*ones((n+1,1))
b[0] = 100
b[n] = 200
T = zeros((n+1,1))

xt = arange(0, l + l/200, l/200)
y = []
for i in xt:
  y.append(-S*i*i+25000*i + 100)
#plot(xt, y, color='red')

T = tdma_solve((A, Ad, Au), b)
plot(x,T)
#Energy balance, left face + right face = S*l
print -(k/dx)*(T[0]-T[1]), -(k/dx)*(T[n]-T[n-1])
print -(k/dx)*(T[0]-T[1]) - (k/dx)*(T[n]-T[n-1])
print S*l
show()
