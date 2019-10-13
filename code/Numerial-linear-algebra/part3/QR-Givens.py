import numpy as np
import copy
import time
start=time.time()

def Given(a,b):
    if b ==0:
        if a>=0:
            c=1
            s=0
        else:
            c=-1
            s=0
    else:
        if abs(b)>abs(a):
            t=a/b
            s=np.sign(b)/np.sqrt(1+t**2)
            c=s*t
        else:
            t=b/a
            c=np.sign(a)/np.sqrt(1+t**2)
            s=c*t
    g=np.array([[c,s],[-s,c]])
    return g

def qr_g(a):
    a=a.astype(float)
    m,n=a.shape
    q=np.eye(m)
    for k in range(n):
        for i in range(m)[k+1:]:
            g=Given(a[k,k],a[i,k])
            a[[k,i],k:]=np.matmul(g,a[[k,i],k:])
            q[:,[k,i]]=np.matmul(q[:,[k,i]],g.T)
    return q,a

a=np.array([[3,5,6],
    [4,8,6],
    [9,2,5],
    [8,6,8]])
q,a=qr_g(a)
print(a,q)
print(np.matmul(q,a))
end=time.time()
print(end-start)
