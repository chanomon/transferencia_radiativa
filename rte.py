#!/usr/bin/env python3

#Copyright (C)  YEAR  YOUR NAME.
#    Permission is granted to copy, distribute and/or modify this document
#    under the terms of the GNU Free Documentation License, Version 1.3
#    or any later version published by the Free Software Foundation;
#    with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
#    A copy of the license is included in the section entitled "GNU
#    Free Documentation License".


from matplotlib import pyplot as plt
import math
from astropy.modeling.blackbody import blackbody_lambda
from tqdm import tqdm
from scipy.interpolate import interp1d
c= 3e10 #cm/s
kB = 1.38e-16 #[eg/K]


#Temperature model [K]
def T(x):
    with open("T.dat","r") as Tf:
        T_dat=Tf.readlines()
    z=[]
    Tr=[]
    
    for line in T_dat:
        if line[0] != '#':
            a,b = line.split()
            z.append(float(a))
            Tr.append(float(b))
    f = interp1d(z,Tr)
    return f(x)
#Density model [p/cm3]
def n(x):
    with open("n.dat","r") as nf:
        n_dat=nf.readlines()
    z=[]
    nr=[]
    for line in n_dat:
        if line[0] != '#':
        	a,b = line.split()
        	z.append(float(a))
        	nr.append(float(b))
    f = interp1d(z,nr)
    return f(x)

#[erg /cm2 sec cm ster]
def S(x,wl):
    return blackbody_lambda(wl*1e8, T(x))*1e8 
#opacity [cm-1]
def k(x,wl):
    
    nu = c/wl
    #Ref Dulk (1985) eq.21
    return 1e5 * 0.2*pow(n(x),2)*pow(T(x),-3/2.)*pow(nu,-2)

#optical depth (adimensional)
def tau(dx,x,wl):
    return (dx/2.)*(k(x-dx,wl)+k(x,wl))

def rayleigh(I,wl):
    return I*pow(wl,4)/(2.0*c*kB)



N = 6.96e2 #Number of points in the raypath
I0= 0.0 #[erg /cm2 sec cm ster]
dx = 1e3#[km]
#wl = 6565.#Amstrongs
nu = 1e8 # [Hz]
wl = c/nu

layers = range(1,int(N+1))
I = I0

X=[]
Y=[]

for i in layers:#tqdm(layers):
    x = float(i) * dx
    I = I*math.exp(-tau(dx,x,wl)) + S(x,wl)*(1-math.exp(-tau(dx,x,wl)))
    X.append(x)
    Y.append(rayleigh(I.value,wl))
    #print(x,I.value)
    #pass
    #print(i,I)
print("%e"%rayleigh(I.value,wl))
fig,ax =plt.subplots()
ax.plot(X,Y)
ax.set_xscale("log")
ax.set_yscale("log")
plt.show()
