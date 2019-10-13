import numpy as np
import copy
import time
start=time.time()

def MGS(A):
    m = len(A)
    n = len(A[0])
    r = np.zeros((n,n))
    q = np.zeros((m,n))
    if any(A[:,0])==0:
        q[:,0]=np.zeros(m)
    else:
        r[0,0]=np.linalg.norm(A[:,0],ord=2)
        q[:,0]=A[:,0]/r[0,0]
    for m in range(n-1):
        k=m+1
        q[:,k]=A[:,k]
        for i in range(k):
            r[i,k]=np.matmul(q[:,i],q[:,k])
            q[:,k]=q[:,k]-r[i,k]*q[:,i]
        if any(q[:,k])!=0:
            r[k,k]=np.linalg.norm(q[:,k],ord=2)
            q[:,k]=q[:,k]/r[k,k]
    return r,q

A=np.array([[1,0,0],
    [0,1,0],
    [0,0,-1],
    [0,-1,0]])
(R,Q)=MGS(A)
print(Q,R)
print(np.matmul(Q,R))
end = time.time()
print(end-start)
