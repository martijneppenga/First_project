import numpy as np
from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython import display
from constants import *


cylinder = True
boundary = compute_structure(N,M,cylinder)
rho = np.ones((N,M))
velocity = 0.001*np.ones((N,M,2))
velocity[:,:,1] = 0 
velocity[boundary,:] = 0


flow_density = equilibrium_function((N,M,9),rho,w,velocity,e)



## Main loop ##
for t in range(Numtimesteps):

    rho = np.sum(flow_density, axis=2)
    velocity = velocity_calc(flow_density,e,rho)
    velocity[:,:,0] +=0.001 * v_int
    f_eq = equilibrium_function((N,M,9),rho,w,velocity,e)
    flow_density_new = relaxation(flow_density,f_eq,tau) 
    flow_density_new = bounce_back(flow_density_new,flow_density,boundary,bounce)
    flow_density = flow_fluid(flow_density_new,e,boundary)
    
    # Make images of velocity in the x and y direction and of the speed
    if (t % 1000 == 0):
        plt.clf()
        fig=plt.figure(figsize=(10, 5), dpi= 80, facecolor='w', edgecolor='k')
        
        fig=plt.figure(figsize=(10, 5), dpi= 80, facecolor='w', edgecolor='k')
        
        
  
        plt.subplot(1, 2, 1)
        plt.imshow((np.sqrt(velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2)).transpose())
        plt.colorbar()
            
        plt.subplot(1, 2, 2)
        plt.plot(velocity[:,50,0])
        plt.tight_layout()
        display.clear_output(wait=True)
        plt.show()
