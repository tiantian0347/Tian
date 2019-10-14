import numpy as np
import time
import copy
start = time.time()
def House(x):
    n = len(x)
    m = 0.
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
        a = np.sqrt(x[0]**2+m)
        if x[0]<0:
            v[0]=x[0]-a
        else:
            v[0]=-(m/(x[0]+a))
        b=2/(v[0]**2+m)
    H = np.eye(n)-(b*np.outer(v,v))
    return b,v

def qr(A):
    A=A.astype(float)
    (m,n) = (len(A),len(A[0]))
    q = np.eye(m)
    r = np.zeros((m,n))
    for k in range(n):
        y = copy.deepcopy(A[k:,k])
        (b,v)=House(y)
        v=np.array([v])
        t=np.matmul(v,A[k:,k:])
        A[k:,k:]-=b*np.matmul(v.T,t)
        q[:,k:]-=b*np.matmul(np.matmul(q[:,k:],v.T),v)
    return q,A

def ort(a):
    k=0
    n,n=a.shape
    p=n-2
    z=np.eye(n,p)
    times=100
    while k<times:
        y=np.matmul(a,z)
        (z,r)=qr(y)
        k+=1
    return z

A=np.array([[5.,4.,3.],
    [0.,4.,5.],
    [1.,0.,-11.]])
p=ort(A)
print(p)
end = time.time()
print(end-start)
