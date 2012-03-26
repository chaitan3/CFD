from pylab import *
from tdma import tdma_solve
n = 10
l = 0.04
dx = l/n
k = 20
h = 10
Tinf = 25
Ar = 2e-5
P = 0.015
#Linear source term
Sc = h*P*Tinf/Ar
Sp = -h*P/Ar
x = arange(0, l + dx, dx)

A =(2*(k/dx)-Sp*dx)*ones((n+1, 1))
A[0] = 1
A[n] = 1
Au = -(k/dx)*ones((n+1, 1))
Au[0] = 0
Ad = -(k/dx)*ones((n+1, 1))
Ad[n] = 0
b = (Sc*dx)*ones((n+1,1))
b[0] = 100
b[n] = 40
T = zeros((n+1,1))

xt = arange(0, l + l/40, l/40)
y = []
m = (h*P/(k*Ar))**0.5
for i in xt:
  t = ((40-Tinf)*sinh(m*i)/(100-Tinf)+sinh(m*(l-i)))/sinh(m*l)
  y.append(Tinf + (100-Tinf)*t)
#plot(xt, y)
T = tdma_solve((A, Ad, Au), b)
plot(x,T)
show()
