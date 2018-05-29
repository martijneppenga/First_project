import numpy as np
import matplotlib.pyplot as plt

from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider
import sys
import matplotlib.patches as mpatches
from matplotlib import cm
from IPython import display
from functions import *
from constants import *


boundary = compute_structure(N,M,structure,r,xcyl,ycyl)
rho = np.ones((N,M))
velocity = np.zeros((N,M,2))


velocity[:,:,0] = init_velocity(maxv,N,M)
velocity[boundary,:] = 0
u = np.copy(velocity)




flow_density = equilibrium_function((N,M,9),rho,w,velocity,e)

fig, ax = plt.subplots(1, figsize=(8,4), facecolor=(1,1,1))
fig.subplots_adjust(left=0, right=1, bottom=0)
xx, yy = np.meshgrid(np.linspace(-2,3,500), np.linspace(-1,2,500))

def make_frame(t):
    global flow_density
    for i in range(100):
        flow_density[-1,:,4] = flow_density[-2,:,4]
        flow_density[-1,:,5] = flow_density[-2,:,5]
        flow_density[-1,:,6] = flow_density[-2,:,6]

        rho = np.sum(flow_density, axis=2)


        velocity = velocity_calc(flow_density,e,rho)
        velocity[0,:,:] = u[0,:,:]
        rho[0,:] = rho_zou_he(flow_density,u)
        f_eq = equilibrium_function((N,M,9),rho,w,velocity,e)
        flow_density = zou_he_boundary(flow_density,f_eq)

        
        
        flow_density_new = relaxation(flow_density,f_eq,tau) 
        flow_density_new = bounce_back(flow_density_new,flow_density,boundary,bounce)



        flow_density = flow_fluid(flow_density_new,e,boundary)
    vabs =  np.sqrt(velocity[:,:,0]**2 + velocity[:,:,1]**2)
    vabs[boundary] = float('NaN')
    
    #vx = rho
    ax.clear()
    ax.axis('off')
    ax.set_title("vabs", fontsize=16)    

    # the varying weights make the points appear one after the other
    ax.imshow(vabs.transpose())

    

    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration = 10)
animation.write_gif("vcylinder_hr.gif", fps=30)