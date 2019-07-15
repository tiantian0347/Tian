import numpy as np
import matplotlib.pyplot as plt
a=0
b=1
ya=1
step=200
y=np.zeros([step+1])
x=np.linspace(a,b,step+1)
y[0]=ya
u=np.exp(x)
for i in range(0,step):
    y[i+1]=y[i]+((b-a)/step)*y[i]
    print(y[i])
plt.xlabel("x")
plt.ylabel("error")
plt.plot(x,y-u) 
plt.show()

