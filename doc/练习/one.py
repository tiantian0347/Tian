"""
五点差分法
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse.linalg import spsolve
#计算系数矩阵
def coefficientmatrix(h1,h2,n):
    n0 = (n-1)*(n+1)
    A = np.zeros((n0,n0),dtype = float)
    c = (2/(h1*h1)+2/(h2*h2))*np.ones(n0)
    d = (-1/(h1*h1))*np.ones(n0-1)
    e = (-1/(h2*h2))*np.ones(n0-n-1)
    A = np.diag(c)+np.diag(d,1)+np.diag(d,-1)+np.diag(e,n+1)+np.diag(e,-n-1)
    for i in range(n-1):
        m = i*(n+1) #(n-1)*(n+1)个元素
        A[m,m] = 1.
        A[m,m+1] = -1.
        A[m,m-1] = 0.#第m+1行除特殊外全为0
        if i>0:
            A[m-n-1,m]=0.
            A[m-n-2,m-1]=0.
            A[m,m-n-1]=0.
            A[m-1,m-n-2]=0.
            A[m-1,m-1] = 1.#对角线
            A[m-1,m-2] = -1.#对角线左边
            A[m-1,m] = 0.
        A[n0-2-n,n0-1]=0
        A[n0-1,n0-2-n]=0
        A[n0-1,n0-1]=1
        A[n0-1,n0-2]=-1
        A[n0-1,n0-n-2]=0
    return A


def fij(n,h1,h2,x,y,x_low,x_up,dy_left,dy_right):
    #方程右端的值包括两部分，一部分为边界上的，另一部分为函数本身
    n0=(n-1)*(n+1)
    f2=np.zeros((n-1,n+1),dtype=float)
    f4=np.zeros((n-1,n+1),dtype=float)
    #计算每个网格上的方程右边函数值f(x,y)
    for i in range(1,n):#沿着y方向，点数少
        for j in range(0,n+1):
            f2[i-1,j]=fxy(x[j],y[i])
            f4[i-1,j]=fab(x[j],y[i])
    f3=f2.flatten()#变成一位数组
    for k in range(n-1):
        m = k*(n+1)
        f3[m] = dy_left*h1
        f3[m-1] = dy_right*h1
#    print(f3)
    f5=f4.flatten()
    return f3,f5#总的函数值

def fxy(x,y):#函数
    f2=np.cos(3*x)*np.sin(np.pi*y)
    return f2

def fab(x,y):
    f1=1/(9+np.pi**2)*np.cos(3*x)*np.sin(np.pi*y)
    return f1


n=4#网格数
a_x=0.#x方向边界
b_x=np.pi
c_y=0.#y方向边界
d_y=1.
h1=(b_x-a_x)/n
h2=(d_y-c_y)/n

A0=coefficientmatrix(h1,h2,n)#计算系数矩阵
#print(A0)


x_cood=np.arange(a_x,b_x+h1,h1)
#n+2个网格，n+2个未知数
y_cood=np.arange(c_y,d_y+h2,h2)
#y有n+1个点，n-1个未知数
x_low=0.
x_up=0.
dy_left=0.
dy_right=0.

(b_matrix,p)=fij(n,h1,h2,x_cood,y_cood,x_low,x_up,dy_left,dy_right)
#计算等式右端矩阵b
f_1d=spsolve(A0,b_matrix)#调用共轭向量法
print('数值解为',f_1d)
print('解析解为',p)

er=np.linalg.norm((f_1d-p).reshape(-1,1),ord=2)/(n**2-1)
print('error为',er)
result=f_1d.reshape(n-1,n+1)#转换成二维矩阵
y_cood2=np.arange(c_y+h2,d_y,h2)#y有n+1个点，n-1个未知数
X, Y = np.meshgrid(x_cood, y_cood2)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, result, cstride=1, cmap=cm.viridis)
plt.show()

