import numpy as np
import matplotlib.pyplot as plt
import datetime
import importlib
import constants; importlib.reload(constants); from constants import *
def start(A,dose_t,dose_spine):
    
    
    #for plot font
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('font', size=16)
    
    #Plot of the 2D projection of the 3D absorption matrix
    plt.figure(figsize=(6,6)) 
    plt.imshow(np.sum(A,axis = 0))
    plt.xlabel('z(cm)')
    plt.ylabel('y(cm)')
    plt.title('Deposited energy vs coordinates')
    plt.colorbar()
    plt.show()
    
    #1D Projection of the 3D absorption matrix
    plt.plot(np.sum(np.sum(A,axis=0),axis = 0))
    plt.xlabel('E (photon weight)')
    plt.ylabel('y(cm)')
    plt.title('Deposited energy vs coordinates')
    plt.show()
    
    #Some properties of the dose deposition in the tissue
    print("Tumour dose = ", dose_t)
    print("Spine dose = ", dose_spine)
    print("Dose tumour/Dose_spine", dose_t/dose_spine)


    
