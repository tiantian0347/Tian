import numpy as np
import copy
import time
start=time.time()

def Given(x,i,j):
    a=copy.deepcopy(x[i-1])
    b=copy.deepcopy(x[j-1])
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
    n=len(x)
    G=np.eye(n)
    G[i-1,i-1]=c
    G[j-1,j-1]=c
    G[i-1,j-i]=s
    G[j-1,i-1]=-s
    return c,s,G


x=np.array([1,2])
(c,s,G)=Given(x,1,2)
y=np.matmul(G,x)
print(G,y)
end=time.time()
print(end-start)
