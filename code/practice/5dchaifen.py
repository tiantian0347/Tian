'''
五点差分
△u=f(x,y),求u在(0,pi)*(0,1)上的数值解
'''
import numpy as np
import copy
import scipy.linalg
def f(x,y):
    z=-np.cos(3*x)*np.sin(np.pi*y)
    return z
#真解
def zj(x,y):
    z=np.cos(3*x)*np.sin(np.pi*y)/(9+np.pi**2)
    return z
#生成三对角矩阵
def sdj(a,b,n):
    x=a*np.identity(n)
    for i in range(n):
        if i!=0:
            x[i,i-1]=b
        if i!=n-1:
            x[i,i+1]=b
    return x
def yyy(n):
    a,b=0,np.pi
    c,d=0,1
    h1=(b-a)/n
    h2=(d-c)/n
    xdate=np.linspace(a,b,n+1)
    ydate=np.linspace(c,d,n+1)
    #边界条件
    ux_0=np.zeros([n+1])
    ux_1=np.zeros([n+1])
    #系数矩阵
    u1=sdj(-1/h1**2-1/h2**2,1/h1**2,n-1)
    u2=sdj(-1/h1**2-1/h2**2,1/h2**2,n-1)
    u0=np.identity(n-1)
    u=np.kron(u0,u2)+np.kron(u1,u0)
    for i in range((n-1)**2):
        if (i+1)%(n-1)<2:
            u[i,i]=-1/h1**2-2/h2**2
    #函数值
    b=np.zeros([n-1,n-1])
    for i in range(n-1):
        for j in range(n-1):
            b[i,j]=f(xdate[i+1],ydate[j+1])
    b[:,0]-=ux_0[1:-1]/h2**2
    b[:,-1]-=ux_1[1:-1]/h2**2
    b=b.flatten()
    #求解方程组
    x=np.linalg.solve(u,b)
    #真解
    kkk=np.zeros([n-1,n-1])
    for i in range(n-1):
        for j in range(n-1):
            kkk[i,j]=zj(xdate[i+1],ydate[j+1])
    kkk=kkk.flatten()
    #误差
    error=np.linalg.norm(x-kkk,ord=2)
    lk=np.matmul(u,x)-b
    x=x.reshape(n-1,n-1)
    b=b.reshape(n-1,n-1)
    return error
er=np.zeros([50])
for i in range(50)[2:]:
    q=yyy(i)
    er[i]=q
print(er)



























