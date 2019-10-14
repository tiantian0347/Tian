import numpy as np
from numpy.linalg import solve,norm
pi=np.pi
sin=np.sin
cos=np.cos
# In[]:right_hand_term
def rhs(x,y):
    f=cos(3*x)*sin(pi*y)
    return f
# In[]:main function
def elliptic(N):
    h1=pi/N
    h2=1/N
    x=np.linspace(0,pi,N+1)
    y=np.linspace(0,1,N+1)
    X,Y=np.meshgrid(x,y)
    B=np.zeros(shape=(N+1,N+1),dtype=float)#
    C=np.zeros(shape=(N+1,N+1),dtype=float)#
    F=(cos(3*X)*sin(pi*Y)).reshape(-1,1)
    for i in range(1,N):
        C[i,i]=-1/(h2**2)
        B[i,i-1]=-1/(h1**2)
        B[i,i]=2/(h2**2)+2/(h1**2)
        B[i,i+1]=-1/(h1**2)
    A=np.kron(np.eye(N+1),B)+np.kron(np.eye(N+1,k=1)+np.eye(N+1,k=-1),C)
# In[]:边界处理
    for i in range(N+1):
        for j in range(N+1):
            if i==j:
                A[i,j]=1
                A[N**2+N+i,N**2+N+j]=1
            else:
                A[i,j]=0
                A[N**2+N+i,N**2+N+j]=0
        F[i]=0
        F[N**2+N+i]=0

    for i in range(N-1):
        A[(i+1)*(N+1),(i+1)*(N+1)]=1
        A[(i+1)*(N+1),(i+1)*(N+1)+1]=-1
        F[(i+1)*(N+1)]=0
    for i in range(N-1):
        A[N+(i+1)*(N+1),N+(i+1)*(N+1)]=1
        A[N+(i+1)*(N+1),N+(i+1)*(N+1)-1]=-1
        F[N+(i+1)*(N+1)]=0
    for i in range(N+1):
        for j in range(N+1):
            A[i,N+j+1]=0
            A[N**2+N+i,N**2+j-1]=0

#系数矩阵右下角分块矩阵化为单位矩阵
    for i in range(N+1):
        for j in range(N+1):
            if i==j:
                A[N**2+N+i,N**2+N+j]=1
            else:
                A[N**2+N+i,N**2+N+j]=0
    uh=solve(A,F)
    uh=np.reshape(uh,(N+1,N+1))
    print(uh)
    exact_u=(1/(9+pi**2))*(cos(3*X)*sin(pi*Y))
    return A,X,Y,uh,exact_u

# In[]:误差
N1=4
A1,X1,Y1,uh1,u_true1=elliptic(N1)
e1=norm((uh1-u_true1).reshape(-1,1),ord=2)/(N1**2)
print('N=%d,l2误差为%.6f'%(N1,e1))
print('数值解为%.6f',uh1)
print('解析解为%.6f',u_true1)
# In[]:
N2=8
A2,X2,Y2,uh2,u_true2=elliptic(N2)
e2=norm((uh2-u_true2).reshape(-1,1),ord=2)/(N2**2)
R1=np.math.log(e1/e2)/np.math.log(2)
print('N=%d,l2误差为%.6f,误差阶为%.6f'%(N2,e2,R1))
# In[]:
N3=16
A3,X3,Y3,uh3,u_true3=elliptic(N3)
e3=norm((uh3-u_true3).reshape(-1,1),ord=2)/(N3**2)
R2=np.math.log(e2/e3)/np.math.log(2)
print('N=%d,l2误差为%.6f,误差阶为%.6f'%(N3,e3,R2))
# In[]:
N4=32  
A4,X4,Y4,uh4,u_true4=elliptic(N4)    
e4=norm((uh4-u_true4).reshape(-1,1),ord=2)/(N4**2)    
R3=np.math.log(e3/e4)/np.math.log(2)                                
print('N=%d,l2误差为%.6f,误差阶为%.6f'%(N4,e4,R3))

# In[]: 
'''                                            
N5=64    
A5,X5,Y5,uh5,u_true5=elliptic(N5)    
e5=norm((uh5-u_true5).reshape(-1,1),ord=2)/(N5**2) 
R4=np.math.log(e4/e5)/np.math.log(2)                                
print('N=%d,l2误差为%.6f,误差阶为%.6f'%(N5,e5,R4)) 
'''
# In[]:Draw picture
'''
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

fig=plt.figure(1)
ax=Axes3D(fig)
A1,X1,Y1,uh1,exact_u1=elliptic(N1)
h1=pi/N1
h2=1/N1
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u(x,y)')
ax.plot_surface(X1,Y1,uh1,rstride=1,cstride=1,camp=plt.get_camp('rainbow'))
plt.title('uh')

fig=plt.figure(2)
ax=Axes3D(fig)
A2,X2,Y2,uh2,exact_u2=elliptic(N2)
h1=pi/N2
h2=1/N2
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u(x,y)')
ax.plot_surface(X2,Y2,uh2,rstride=1,cstride=1,camp=plt.get_camp('rainbow'))
plt.title('uh')

fig=plt.figure(3)
ax=Axes3D(fig)
A3,A3,Y3,uh3,exact_u3=elliptic(N3)
h1=pi/N3
h2=1/N3
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u(x,y)')
ax.plot_surface(X3,Y3,uh3,rstride=1,cstride=1,camp=plt.get_camp('rainbow'))
plt.title('uh')

fig=plt.figure(4)
ax=Axes3D(fig)
A4,X4,Y4,uh4,exact_u4=elliptic(N4)
h1=pi/N4
h2=1/N4
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u(x,y)')
ax.plot_surface(X4,Y4,uh4,rstride=1,cstride=1,camp=plt.get_camp('rainbow'))
plt.title('uh')

fig=plt.figure(5)
ax=Axes3D(fig)
A5,X5,Y5,uh5,exact_u5=elliptic(N5)
h1=pi/N5
h2=1/N5
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('u(x,y)')
ax.plot_surface(X5,Y5,exact_u5,rstride=1,cstride=1,camp=plt.get_camp('rainbow'))
plt.title('exact_u5')

plt.show()
'''
