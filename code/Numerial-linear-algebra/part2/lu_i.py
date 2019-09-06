import numpy as np
def bd_lu(a):
    (n,m)=np.shape(a)
    for k in range(0,n):
        for i in range(k+1,n):
            a[i,k]=a[i,k]/a[k,k]
            for j in range(k+1,n):
                a[i,j]=a[i,j]-a[i,k]*a[k,j]
    return a
def main():
    A=np.array([[2.,2.,3.],
        [4.,7.,7.],
        [-2.,4.,5.]])
    print(bd_lu(A))

main()
