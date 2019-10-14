import nump as np

def par(a,h,m,x_low,x_up,t_left,t_right):
    J=(x_up-x_low)/h
    N=(t_right-t_left)/m
    n0=(J-1)(N-1)
    r=a*m/(h**2)
    c = (1+2*r)*np.ones(n0)
    d = (-r)*np.ones(n0-1)
    A = np.diag(c)+np.diag(d,1)+np.diag(d,-1)
    for i in range(0,N-1):
        j=i*(J-1)
        A[j,j-1] = 0
        A[j-1,j] = 0
        if i>0:
            A[j,j-J] = -1
    return A
def f(


