import numpy as np
import matplotlib.pyplot as plt
import time
import numba
import sys
import importlib
import functions; importlib.reload(functions); from functions import *
import constants; importlib.reload(constants); from constants import *
def start():
    #Initial conditions
    w_refl = 0
    r0,u0 = initialization(nbps, N, L, a, mode)
    #initialisation  
    w_refl = 0
    w = np.ones((N,))

    A = np.zeros((Nbins,Nbins,Nbins))                       #Matrix which tracks where the dose was deposited
    T = np.zeros((Nbins,Nbins,Nbins))                       #Fractional transport, probably matrix with same size as our grid


    mu_t = mu_a + mu_s  #attenuation coefficient

    r = np.copy(r0)
    u = np.copy(u0)



    if abs(g)>1:
        sys.exit("g (anisotropy of scattering) is too large, please choose g such that |g|<1")
    if 1/mu_t*4>L:
        sys.exit("The mean free path lenght is too large(1/mu_t). Please choose a smaller mean free"+
                 " path length or large box size (L)")    

    start = time.time()   
    while np.sum(w)>0:




        r_old = np.copy(r)

        #hop
        r,s = hop(r,u,mu_t)
        if cavity:       
            #Reflect and transmit photons that travel through the cavity
            r,u,w,s,out = cavity_trans_refl(r,r_old,u,w,s,cavcent,cavr,n1,n2,mu_t)
        else:
            out = np.ones(len(w),dtype=bool)
        #reflection on boundary of the system
        r[:,out],u[:,out],w[out],w_refl = reflection(r[:,out],u[:,out],s[out],L,w[out],w_refl,n1,n2)

        if np.sum(out)>0:
            #absorption
            A,w[out] = absorb(A,r[:,out],Nbins,dL,w[out],mu_a,mu_t)

            #scattering
            u[:,out] = scattering(u[:,out],g)

        #terminate?
        w,u,r = terminate(w,u,r,threshold,p_term)
    x = np.linspace(0,Nbins-1,Nbins)
    y = np.copy(x)
    z = np.copy(x)
    [X,Y,Z] = np.meshgrid(x,y,z)
    dose_spine = np.sum(A[((X-spine_c[0]/dL)**2 + (Z-spine_c[1]/dL)**2 <spine_r**2)])/(4/3*np.pi*t_r**3)                #Dose in the spine
    dose_t = np.sum(A[((X-t_c[0]/dL)**2 + (Y-t_c[1]/dL)**2 + (Z-t_c[2]/dL)**2)<t_r**2])/(20*np.pi*spine_r**2)           #Dose in the tumor

    end = time.time()
    dt = end-start
    print("Time elapsed", dt)

    return A,dose_t,dose_spine




