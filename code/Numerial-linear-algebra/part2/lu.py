import numpy as np


def bd_LU(a):
    (n,m)=np.shape(a)
    l=np.identity(n)
    u=np.zeros(shape=(n,n))

    for k in range(0,n):
        for i in range(k+1,n):
            l[i,k]=a[i,k]/a[k,k]
        for j in range(k,n):
            u[k,j]=a[k,j]
        for i in range(k+1,n):
            for j in range(k+1,n):
                a[i,j]=a[i,j]-l[i,k]*u[k,j]
    return l,u

def for_tri(b,l):
    (n,m)=np.shape(l)
    y=np.zeros(n)
    y[0]=b[0]
    for i in range(1,n):
        for j in range(0,i):
            b[i]=b[i]-l[i,j]*y[j]
        y[i]=b[i]
    return y

def back_tri(y,u):
    (n,m)=np.shape(u)
    x=np.zeros(n)
    x[n-1]=y[n-1]/u[n-1,n-1]
    for i in range(0,n-1):
        i=n-2-i
        for j in range(0,i+1):
            y[j]=y[j]-u[j,i+1]*x[i+1]
        x[i]=y[i]/u[i,i]
    return x
def main():
    A=np.array([[2.,2.,3.],
        [4.,7.,7.],
        [-2.,4.,5.]])
    c=np.array([3,3,2])
    l,u=bd_LU(A)
    print(l)
    print(u)
    y=for_tri(c,l)
    x=back_tri(y,u)
    print(y)
    print(x)
    print(np.dot(l,y))
    print(np.dot(u,x))
    print(np.dot(l,np.dot(u,x)))

main()
