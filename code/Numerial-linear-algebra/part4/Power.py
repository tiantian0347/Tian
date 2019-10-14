'幂迭代'

import numpy as np
import copy

def power(a):
    k=0
    n,n=a.shape
    x=np.ones(n)/np.sqrt(n)
    y1,y2=0,1
    while abs(y2-y1)>error:
        y1=copy.deepcopy(y2)
        t=np.matmul(a,x)
        x=t/np.linalg.norm(t,ord=2)
        y2=np.matmul(x,np.matmul(a,x))
        k+=1
    return y2,x

a=np.array([[1,0,0],
    [0,3,0],
    [0,0,5]])
error=0.000000000001
c,x=power(a)
print(c,x)
