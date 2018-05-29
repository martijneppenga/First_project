import numpy as np

# Flow constants
#---------------
Numtimesteps = 100000        # amount of cycles

N = 400                            # Lattice points in x-direction
M = 200                            # Lattice points in y-direction
maxv = 0.01                        #maximum velocity of the laminar inflow
nu       = 0.01
tau      = (6*nu + 1)/2
#v_int = 0.01                      # maximum velocity of Poiseuille flow

#Structure inside the pipe
#---------
structure = "maze"                 #choose from "maze", "cylinder", "none"

xcyl = N/4
ycyl = M/2
r = 5

w        = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36])                         #weights for the different directions
e        = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]])     #directions in which the lquid can flow
bounce = np.array([[0,1,2,3,4,5,6,7,8],[0,5,6,7,8,1,2,3,4]])                           #Array that pairs the opposite directions of e
