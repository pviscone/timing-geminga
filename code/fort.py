import numpy as np
import ztest as zt
ph_time=np.load("ph_time.npy")/1000000

#ph_time=np.linspace(-10,10,11)
f0_start=4.2175
f1_start=-1.95250e-1

n0=50
delta0=2

n1=50
delta1=0.01

f0_arr=np.linspace(f0_start-delta0,f0_start+delta0,n0)
f1_arr=np.linspace(f1_start-delta1,f1_start+delta1,n1)

d=zt.mod.matrix(ph_time,f0_arr,f1_arr)
print(d)