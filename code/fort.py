""" import numpy as np
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
print(d) """

#arr1= array,matrice fasi di ogni fotrone
import numpy as np

def arr(f0,f1):
    a=np.array([1,2,3])
    return  f0*a+f1


def arr_funct(arr1,arr2=np.arange(4,6)):
    arr2=np.reshape(arr2,(arr2.shape[0],1,1))
    return arr2*arr1


def wrapper(f0_arr,f1_arr):
    return  arr_funct([[arr(f0,f1) for f0 in f0_arr] for f1 in f1_arr])


""" mat=np.array([[1,3,5],[2,4,1],[1,1,1],[1,1,1]])
cc=arr_funct(mat) """

a=np.arange(10,15)
b=np.arange(100,109)


wrapper(a,b)