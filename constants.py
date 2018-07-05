import numpy as np
#Constants
#---------

#simulation constants
Nbins = 100                              #number of bins
L = 100                                  #box size
threshold = 1E-4                         #threshold value of the weights to terminate photon
p_term = 0.1                             #probability of terminating photon

#material constants
mu_a = 0.03                              #absorption coefficient
mu_s = 0.04                              #scattering coefficient
g = 0.9                                  #anisotropy of scattering
n1 = 1.5                                 #refractive index inside medium 
n2 = 1                                   #refractive index outside medium

#beam properties
N = int(1e4)                             #number of photons
a = 10                                   #beam radius (circular)
P = 1                                    #beam power
mode = "single"                          #beams from one side "single" or multiple side "multiple"
nbps = 40                                #number of beams per side

#cavity
cavity = False                           #cavity in simulation
cavcent = np.array([[L/2],[L/2],[L/4]])  #coordinates of cavity
cavr = L/5                               #radius cavity

#Tissue properties
dL = L/Nbins                             #Compute length scale of the bins in A
t_c = np.array([L/2,L/2,L/2])/dL         #Tumor location
t_r = 2.5/dL                             #Tumor radius (sphere)
spine_c = np.array([L/3,L/2])/dL         #Spine location
spine_r = 1/dL                           #Spine radius (cylinder)
