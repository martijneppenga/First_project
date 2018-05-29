import numpy as np
from functions import *
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython import display
from constants import *
import time

boundary = compute_structure(N,M,structure,r,xcyl,ycyl)
rho = np.ones((N,M))
velocity = np.zeros((N,M,2))

velocity[:,:,0] = init_velocity(maxv,N,M)
velocity[boundary,:] = 0
u = np.copy(velocity)


flow_density = equilibrium_function((N,M,9),rho,w,velocity,e)



## Main loop ##
for t in range(Numtimesteps):
    flow_density[-1,:,1] = flow_density[-2,:,1]
    flow_density[-1,:,2] = flow_density[-2,:,2]
    flow_density[-1,:,8] = flow_density[-2,:,8]
    
    rho = np.sum(flow_density, axis=2)
 
    
    
    
    
    velocity = velocity_calc(flow_density,e,rho)
    velocity[0,:,:] = u[0,:,:]
    
       
    #velocity[0,:,:] = u[0,:,:] 
    #velocity[:,:,0] +=0.001 * v_int
    f_eq = equilibrium_function((N,M,9),rho,w,velocity,e)
    flow_density_new = relaxation(flow_density,f_eq,tau) 
    flow_density_new = bounce_back(flow_density_new,flow_density,boundary,bounce)
    


    
    
    rho_temp = flow_density[0,:,0]+flow_density[0,:,3]+flow_density[0,:,7]+2*(flow_density[0,:,1]+flow_density[0,:,2]+flow_density[0,:,8])
    flow_density[0,:,5] = flow_density[0,:,1]-2/3*rho_temp*u[0,:,0]
    flow_density[0,:,6] = flow_density[0,:,2]+0.5*(flow_density[0,:,3]-flow_density[0,:,7])-1/2*rho_temp*u[0,:,0]
    flow_density[0,:,4] = flow_density[0,:,8]-0.5*(flow_density[0,:,3]-flow_density[0,:,7])-1/2*rho_temp*u[0,:,0]

    
    flow_density = flow_fluid(flow_density_new,e,boundary)

    
    #flow_density[-1,:,6] = 0
    #flow_density[-1,:,5] = 0
    #flow_density[-1,:,4] = 0

    
    # Make images of velocity in the x and y direction and of the speed
    if (t % 100 == 0):
        plt.clf()
        fig=plt.figure(figsize=(10, 5), dpi= 80, facecolor='w', edgecolor='k')
        
        fig=plt.figure(figsize=(10, 5), dpi= 80, facecolor='w', edgecolor='k')
        
        
  
        plt.subplot(1, 2, 1)
        plt.imshow((np.sqrt(velocity[:,:,0] ** 2 + velocity[:,:,1] ** 2)).transpose())
        plt.colorbar()
            
        plt.subplot(1, 2, 2)
        #plt.plot(velocity[1,:,0])
        plt.imshow(rho)
        plt.colorbar()
        plt.tight_layout()
        
        display.clear_output(wait=True)
        plt.show()
        time.sleep(1)
