#!/usr/bin/env python3

#Copyright (C)  YEAR  YOUR NAME.
#    Permission is granted to copy, distribute and/or modify this document
#    under the terms of the GNU Free Documentation License, Version 1.3
#    or any later version published by the Free Software Foundation;
#    with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
#    A copy of the license is included in the section entitled "GNU
#    Free Documentation License".



import math
from astropy.modeling.blackbody import blackbody_lambda
from tqdm import tqdm
c= 3e10 #cm/s
kB = 1.38e-16 #[eg/K]


#Temperature model [K]
def T(x):
    return 1e6

#Density model [p/cm3]
def n(x):
    return 1e7


#[erg /cm2 sec cm ster]
def S(x,wl):
    return blackbody_lambda(wl, T(x))
#opacity [cm-1]
def k(x,wl):
    
    nu = c/wl
    #Ref Dulk (1985) eq.21
    return 0.2*pow(n(x),2)*pow(T(x),-3/2.)*pow(nu,-2)

#optical depth (adimensional)
def tau(dx,x,wl):
    return (dx/2.)*(k(x-dx,wl)+k(x,wl))

def rayleigh(I,wl):
    return I*pow(wl,4)/(2.0*c*kB)



N = 6.96e3
I0= 0.0 #[erg /cm2 sec cm ster]
dx = 100e5#[cm]
#wl = 6565.#Amstrongs
nu = 1e8 # [Hz]
wl = c/nu

layers = range(1,int(N+1))
I = I0
for i in tqdm(layers):
    x = float(i) * dx
    I = I*math.exp(-tau(dx,x,wl)) + S(x,wl)*(1-math.exp(-tau(dx,x,wl)))
    pass
    #print(i,I)
print("%e"%rayleigh(I.value,wl))
