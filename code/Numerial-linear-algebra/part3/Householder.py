import numpy as np
import time
import copy
start = time.time()
def House(x):
    n = len(x)
    m = 0
    for i in range(n-1):
        m = x[i+1]**2+m
    v = np.zeros(n,dtype=float)
    v = copy.deepcopy(x)
    if m == 0:
        if x[0]<0:
            v[0] = 2*x[0]
            b = 2/v[0]**2
        else:
            v[0]=0
            b = 0
    else:
        a = np.sign(x[0])*np.sqrt(x[0]**2+m)
        if x[0]<0:
            v[0]=x[0]-a
        else:
            v[0]=-(m/(x[0]+a))
        b=2/(v[0]**2+m)
    H = np.eye(n)-(b*np.outer(v,v))
    return H

x = np.array([5.,7.,3.])
A = House(x)
print(A)
c = np.matmul(A,x)
print(c)
end = time.time()
print(end-start)
