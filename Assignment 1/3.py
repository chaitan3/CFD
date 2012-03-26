from pylab import *
from tdma import tdma_solve
n = 50
#The 2 domains with different conductivity
la = 0.4
lb = 0.2
l = la+lb
dx = l/n
ka = 10
kb = 0.1
h_r = 10.0
h_l = 25.0
T_l = 800.0
T_r = 20.0
x = arange(0, l + dx, dx)

#Index of x-coordinate array that is just before l=0.4
i = int(la*n/l)
#Interface conductivity, harmonic mean
k = 2*ka*kb/(ka+kb)
A = 2*ones((n+1, 1))
A[0] = ka/dx + h_l
#Interface linear equation with different conductivities on left cell and interface
A[i] = k+ka
#Interface linear equation with different conductivities on right cell and interface
A[i+1] = k+kb
A[n] = kb/dx + h_r
Au = -1*ones((n+1, 1))
Au[0] = -ka/dx
Au[i] = -k
Au[i+1] = -kb
Ad = -1*ones((n+1, 1))
Ad[n] = -kb/dx
Ad[i] = -ka
Ad[i+1] = -k
b = zeros((n+1,1))
b[0] = h_l*T_l
b[n] = h_r*T_r
T = zeros((n+1,1))

b1 = 785.688 
A1 = (25*b1-20000)/10
A2 = 100*A1
b2 = b1-39.6*A1
xt = arange(0, l + l/200, l/200)
y = []
for j in xt:
  if j < 0.4:
    y.append(A1*j+b1)
  else:
    y.append(A2*j+b2)
#plot(xt, y)
T = tdma_solve((A, Ad, Au), b)
plot(x,T)
k1 = ka/(0.4-x[i])
k2 = kb/(x[i+1]-0.4)
#Interface Temperature calculated using the specified formula
print (T[i]*k1 + k2*T[i+1])/(k1+k2)
print A1*0.4+b1
show()
