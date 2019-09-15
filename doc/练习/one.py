"""
五点差分法
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import*
#计算系数矩阵
def coefficientmatrix(h1,h2,n):
    n0 = (n-1)*(n+2)
    A = np.zeros((n0,n0),dtype = float)
    i = np.arange(0,n0)#设置中央对角线
    A[i,i]=2/h1+2/h2
    j = np.arange(0,n0-1)#设置里面的两条斜对角线
    A[j,j+1] = -1/h1
    A[j+1,j] = -1/h1
    k = np.arange(0,n0-n-2)#设置外面的两条斜对角线
    A[k,k+n+2] = -1/h2
    A[k+n+2,k] = -1/h2
    for mm in range(0,n-1):#左边界
        m=mm*(n+2)#(n-1)*(n+2)个元素
        A[m,m]=1./h1
        A[m,m-1]=0.
        if m+n+2<n0:
            A[m,m+n+2]=0.
        else:
            A[m,m-n-2]=0.
        if mm>0:
            A[m,m-n-2]=0.
        if mm>0:#右边界
            A[m-1,m-1]=1./h1#对角线
            A[m-1,m-2]=-1./h1#对角线左边
            A[m-1,m]=0.#对角线右边
        if m+n+1<n0:
            A[m-1,m+n+1]=0.#这一行右边次对角线元素为0
            A[m-1,m-n-3]=0.
        else:
            A[m-1,m-n-3]=0.#这一行左边次对角线元素为0
        A[n0-1,n0-1]=1/h1
        A[n0-1,n0-2]=-1/h1
        A[n0-1,n0-n-3]=0
        return A
def fij(n,h1,h2,x,y,x_low,x_up,dy_left,dy_right):
    #方程右端的值包括两部分，一部分为边界上的，另一部分为函数本身
    n0=(n-1)*(n+2)
    f=np.zeros(n0,dtype=float)
    f[-1]=dy_right#最后一个点 下面机航线设置边界上为0的点
    m=np.arange(0,n0,n+2)
    f[m]=dy_left#左边界
    m2=np.arange(n+1,n0,n+2)#右边界
    f[m2]=dy_right
    m3=np.arange(1,n+1)#下边界
    f[m3]=x_low
    m4=np.arange(n0-n-1,n0-1)
    f[m4]=x_up
    f2=np.zeros((n-1,n+2),dtype=float)
    #计算每个网格上的方程右边函数值f(x,y)
    for i in range(1,n):#沿着y方向，点数少
        for j in range(0,n+2):
            f2[i-1,j]=fxy(x[j],y[i])
    f3=f2.flatten()#变成一位数组
    return f3+(f/h1)/h2#总的函数值
def fxy(x,y):#函数
    f2=-math.cos(3*x)*math.sin(pi*y)
    return f2
def cg(A,b,x):#共轭斜量法
    r=b-np.dot(A,x)#r=b-Ax  r也是梯度方向
    p=np.copy(r)
    i=0
    while(max(abs(r))>1.e-10 and i<100):
        print('i',i)
        print('r',r)
        pap=np.inner(np.dot(A,p),p)
        if pap==0:
            return x
        print('pap=',pap)
        alpha=np.inner(r,r)/pap#直接套用公式
        x1=x+alpha*p
        r1=r-alpha*np.dot(A,p)
        beta=np.inner(r1,r1)/np.inner(r,r)
        p1=r1+beta*p
        r=r1
        x=x1
        p=p1
        i=i+1
    return x
n=16#网格数
a_x=0.#x方向边界
b_x=pi
c_y=0.#y方向边界
d_y=1.
h1=(b_x-a_x)/n
h2=(d_y-c_y)/n

A0=coefficientmatrix(h1,h2,n)#计算系数矩阵
print(A0)

uij=np.zeros((n-1)*(n+2),dtype=float)#未知数，初始化


x_cood=np.arange(a_x-h1/2.,b_x+h1,h1)
#n+2个网格，n+2个未知数
y_cood=np.arange(c_y,d_y+h2,h2)
#y有n+1个点，n-1个未知数
x_low=0.
x_up=0.
dy_left=0.
dy_right=0.

b_matrix=fij(n,h1,h2,x_cood,y_cood,x_low,x_up,dy_left,dy_right)
#计算等式右端矩阵b
f_1d=cg(A0,b_matrix,uij)#调用共轭向量法
result=f_1d.reshape(n-1,n+2)#转换成二维矩阵
plt.figure()
y_cood2=np.arange(c_y+h2,d_y,h2)#y有n+1个点，n-1个未知数
print(y_cood2.shape)
contourf(x_cood,y_cood2,result,80,cmap='seismic')
plt.colorbar()
plt.show()
plt.close()
