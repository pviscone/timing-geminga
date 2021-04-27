import numpy as np
import zt
a=np.linspace(250991900.,250994001.,30957)
b=np.linspace(4,5,500)
c=np.linspace(1e-12,1e-13,50)
d=zt.mod.matrix(a.T,b.T,c.T)
print(d)