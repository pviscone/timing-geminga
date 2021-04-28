import os 
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import ztest
import importlib as il

src_dir = os.getcwd()
main_dir = os.path.dirname(src_dir)
data_dir=os.path.join(main_dir, "data")
results_dir = os.path.join(main_dir, "results")

f0_start=4.2175
f1_start=-1.95250e-13
t0= 250992001

class Data():
    def __init__(self,fit_name):
        self.fit_file=fits.open(os.path.join(data_dir,fit_name))
        self.data=Table(self.fit_file[1].data)
        self.gti=Table(self.fit_file[2].data)
        self.data_header=self.fit_file[1].header
        self.gti_header=self.fit_file[2].header
        self.data.add_column(0,name="phase")
        
    def __repr__(self):
        return str(self.data)
    
    def mapping(self,energy,vmax=150):
        gamma_E=self.data[self.data["ENERGY"]>energy]
        plt.figure(figsize=(6,10))
        plt.subplot(211)
        plt.hexbin(gamma_E["L"],gamma_E["B"],cmap='plasma',vmax=vmax)
        plt.colorbar(label='counts')
        plt.xlabel("L [°]")
        plt.ylabel("B [°]")
        plt.title(f"Gamma with Energy>{energy} map")
        
        plt.subplot(212,projection='aitoff')
        l=np.array([x if x<180 else x-360 for x in gamma_E[gamma_E["L"]>180]["L"]])*np.pi/180
        plt.plot(l[0],gamma_E["B"][0]*np.pi/180,".",markersize=20)
        #plt.xlim(-180,180)
        #plt.ylim(-90,90)
        plt.grid()
        plt.show()
    



tt=Data("PSRJ_Geminga_3deg_100mev_tt.fits")
bary=Data("PSRJ_Geminga_3deg_100mev_bary.fits")


tt.mapping(300)
tt.mapping(1000)



plt.figure()
hist=plt.hist(tt.data["ENERGY"],bins=200,range=(0,8000))
plt.figure()
def exp(x,a,t):
    return a*np.exp(-x/t)
popt,pcov=curve_fit(exp,hist[1][3:-1],hist[0][3:],p0=[2000,1000])
x=np.linspace(60,8000,400)
plt.plot(x,exp(x,*popt))
plt.plot(hist[1][3:-1],hist[0][3:],".")



#il.reload(ztest)
n0=500
delta0=0.2

n1=50
delta1=5e-15
""" 
f0_start=2.4179008016032064
f1_start=4.374030612244898e-13
"""

"""
f0_start=2.7595841683366733
f1_start=4.337295918367347e-13
"""

f0_start=2.856978957915832
f1_start=4.2872959183673473e-13

f0_arr=np.linspace(f0_start-delta0,f0_start+delta0,n0)
f1_arr=np.linspace(f1_start-delta1,f1_start+delta1,n1)

ph_time=np.array(bary.data["TIME"])-t0
np.save("ph_time.npy",ph_time)


z_test=ztest.mod.matrix(ph_time,f0_arr,f1_arr)
z_test=(2/len(bary.data))*z_test
#np.save("z_test.npy",z_test)

plt.figure(figsize=(10,10))

#ATTENZIONE AGLI ASSI FORSE SBAGLIATI
plt.imshow(z_test,aspect='auto',extent=[f1_arr[0],f1_arr[-1],f0_arr[-1],f0_arr[0]])
plt.colorbar()
plt.xlabel("$f_1$ [Hz]")
plt.ylabel("$f_0$ [Hz]")
plt.title("Z Test($f_0$ , $f_1$)")

z_test_max=np.max(z_test)
arg=np.argmax(z_test)
f_best_index=np.unravel_index(arg,z_test.shape)
f0_best=f0_arr[f_best_index[0]]
f1_best=f1_arr[f_best_index[1]]









""" #f0 scorre sulle righe ->, f1 scorre sulle colonne |
ffdot=np.array([[(x,y) for x in f0_arr] for y in f1_arr])

def phi(t,f0,f1):
    return (f0*(t-t0)+(f1*(t-t0)**2)/2)%1

def z(k,phi_arr):
    return np.sum(np.cos(k*phi_arr*2*np.pi)**2)+np.sum(np.sin(k*phi_arr*2*np.pi)**2)

def ztest(phi_arr,n=10,N=len(bary.data)):
    return (2/N)*np.sum((z(k,phi_arr) for k in range(1,n+1)))


#ztest=[[[phi(t,f0,f1) for f0 in f0_arr]for f1 in f1_arr] for t in bary.data["TIME"]]


z_mat=[[ztest(np.array([phi(t,x,y) for t in bary.data["TIME"]])) for x in f0_arr] for y in f1_arr]


"""




""" 
np.save("phase.npy",phase)
#phase=np.load("phase.npy",allow_pickle=True) """