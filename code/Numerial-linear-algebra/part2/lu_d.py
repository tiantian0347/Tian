import numpy as np
def d_lu(a):
    (n,m)=np.shape(a)
    for k in range(0,n):
        for j in range(k,n):
            for i in range(0,k):
                a[k,j]=a[k,j]-a[k,i]*a[i,j]
        for i in range(k+1,n):
            for j in range(0,k):
                a[i,k]=a[i,k]-a[i,j]*a[j,k]
            a[i,k]=a[i,k]/a[k,k]
    return a

def main():
    A=np.array([[2,2,3],
        [4,7,7],
        [-2,4,5]])
    print(d_lu(A))

main()
