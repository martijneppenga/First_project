import numpy as np
import matplotlib.pyplot as plt


def initialization(nbps, N, L, a, mode):
    """initialization(nbps, N, L, a, mode):
    Description:
    ------------
    This function creates the initial setup of the incoming photon beams by evenly distributing the beams over the sides of the structure
    
    Parameters:
    -----------
    nbps: int
        Integer indicating the number of beams per side
    N   : int
        Integer indicating the total number of photons
    L   : float
        Length of one side
    a   : float
        Beam radius
    mode : string
        String which defines if a single beam, or multiple beams are used (input: "single", "multiple")
   
    Results:
    --------
     u0: array of size (3,N)
        array containing the initial unit velocity vector for all the N photons with the rows containing the x,y,z-velocity component.
     r0: array of size (3,N)
        array containing the initial positions for all the N photons with the rows containing the x,y,z-velocity component.
    
    """
    #Creating random radius and phi for a random perturbation (based on the beam radius)
    RND = np.random.rand(N,)
    r1 = a*np.sqrt(RND)
    RND = np.random.rand(N,)
    phi_0 = RND*2*np.pi
    
    if mode == "single":
        r0 = 1E-5*np.ones((3,N))                     #3xN matrix of the initial position for every particle
        r0[0,:] = L/2 + r1*np.cos(phi_0)
        r0[1,:] = L/2 + r1*np.sin(phi_0)

        u0 = np.zeros((3,N))   
        u0[2,:] = 1

    if mode == "multiple":
        
    
    
        pbm  = int(N/(4*nbps))                    #photons per beam
        r0 = np.zeros((3,N))                      
        perturb = np.zeros((3,N))                 #initializing perturbation matrix
        r0[0,:] = L/2                             #x-component is always L/2

        center = L/2*np.ones((3,N))
        for i in range(nbps):
            #Create the positions of the beam centers of the first two sides
            r0[1,i*pbm:(i+1)*pbm] = (i+1)/(nbps + 1)*L 
            r0[2,nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm] = (i+1)/(nbps + 1)*L


            #Create the pertubation to account for a finite beam radius
            perturb[1,i*pbm:(i+1)*pbm] = r1[i*pbm:(i+1)*pbm]*np.cos(phi_0[i*pbm:(i+1)*pbm])
            perturb[0,i*pbm:(i+1)*pbm] = r1[i*pbm:(i+1)*pbm]*np.sin(phi_0[i*pbm:(i+1)*pbm])
            perturb[2,nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm] = r1[nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm]*np.cos(phi_0[nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm])
            perturb[0,nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm] = r1[nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm]*np.sin(phi_0[nbps*pbm+i*pbm:nbps*pbm+(i+1)*pbm])        
        for i in range(nbps):
            #Create the positions of the beam centers of the other two sides
            r0[1,2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm] = (i+1)/(nbps + 1)*L
            r0[2,3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm] = (i+1)/(nbps + 1)*L
            r0[1,3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm] = L
            r0[2,2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm] = L

            #Create the pertubation to account for a finite beam radius
            perturb[1,2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm] = r1[2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm]*np.cos(phi_0[2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm])
            perturb[0,2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm] = r1[2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm]*np.sin(phi_0[2*nbps*pbm+i*pbm:2*nbps*pbm+(i+1)*pbm])
            perturb[2,3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm] = r1[3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm]*np.cos(phi_0[3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm])
            perturb[0,3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm] = r1[3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm]*np.cos(phi_0[3*nbps*pbm+i*pbm:3*nbps*pbm+(i+1)*pbm])

        u0 = center-r0                               #Initial velocity points to center
        lengthu0 = np.sqrt(np.sum(u0**2,axis = 0))   #Computing the length
        u0 = u0/lengthu0                             #Normalizing u0
        r0 = r0 + perturb                            #Add the perturbations to account for finite beam radius

    #u0 = np.zeros((3,N))
    #u0[2,:] = 1
    #r0 = np.zeros((3,N))
    #r0[0,:] = L/2 + r1*np.cos(phi_0)
    #r0[1,:] = L/2 + r1*np.sin(phi_0)
    
    return r0,u0

def scattering(u,g):   
    """scattering(u,g):
    Description:
    ------------
    This function scatters the photons in a random direction. It uses the Henyey-Greenstein function to calculate the
    elevation angle scattering. The azimuth scattering angle is uniform 

    Parameters:
    -----------
    u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity component
    g: float
        anisotropy of scattering (needs to be between -1 and 1)

    Results:
    --------
     u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity component
        after the scattering event

    """ 
    #Henyey-Greenstein function for elevation angle scattering for sin(theta) and cos(theta)
    N = np.shape(u)
    N = N[1]
    RND = np.random.rand(N,)

    if g == 0:
        costheta = 2*RND-1
    else:
        t1 = (1-g**2)/(1-g+2*g*RND)
        costheta = (1+g**2-t1**2)/(2*g)

    sintheta = np.sqrt(1-costheta**2)

    #azimuth angle (phi) between scattered and incoming velocity
    RND = np.random.rand(N,)
    phi = RND*2*np.pi                       
    cosphi = np.cos(phi)
    sinphi = np.sin(phi)

    uxx = np.zeros(N)
    uyy = np.copy(uxx)
    uzz = np.copy(uxx)

    #avoid dividing by zero
    div_zero_true  = abs(u[2,:]) >= (1-1e-12)   #divide by zero locations in live array
    div_zero_false = div_zero_true == False       
    lx = div_zero_false #just a short notation such that the long name div_zero_false is not always neccasary to use


    #Calculate new trajectory vector elements  
    tmp = np.sqrt(1-u[2,div_zero_false]**2)
    uxx[div_zero_false] = sintheta[lx]*(u[0,lx]*u[2,lx]*cosphi[lx]-u[1,lx]*sinphi[lx])/tmp + u[0,lx]*costheta[lx]
    uyy[div_zero_false] = sintheta[lx]*(u[1,lx]*u[2,lx]*cosphi[lx]+u[0,lx]*sinphi[lx])/tmp + u[1,lx]*costheta[lx]
    uzz[div_zero_false] = -sintheta[lx]*cosphi[lx]*tmp + u[2,lx]*costheta[lx]

    uxx[div_zero_true] = sintheta[div_zero_true]*cosphi[div_zero_true]
    uyy[div_zero_true] = sintheta[div_zero_true]*sinphi[div_zero_true]
    uzz[div_zero_true] = costheta[div_zero_true]*np.sign(u[2,div_zero_true])

    u = np.array([uxx, uyy, uzz])
    return u

def absorb(A,r,Nbins,dL,w,mu_a,mu_t):
    """absorb(A,r,Nbins,dL,w,mu_a,mu_t):
    Description:
    ------------
    This function removes a fraction mu_a/mu_t of the photons and puts it in their respective bins (corresponding to their location).
    This fraction is then removed from the photon weights. If the photon is absorbed at an index which does not 1:1 correspond to a bin,
    it is placed in the closest bin instead.

    Parameters:
    -----------
    A: array of size (Nbins,Nbins,Nbins)
       array containing information about where the energy from the photons is absorbed before the called absorption event.
    r: array of size (3,N)
       array with all the x,y,z-locations of each photon
    Nbins: int
       Number of bins. The absorption space is divided into Nbins parts.
    dL: float
       Spacing between bins in units of length.
    w: array of size (N,)
       Contains the weights of the photons that are still propagating, before the absorption event.
    mu_a: float
       Absorption coefficient in units of 1/length.
    mu_t: float
       Total interaction coefficient in units of 1/length


    Results:
    --------
    A: array of size (Nbins,Nbins,Nbins)
       array containing information about where the energy from the photons is absorbed after the called absorption event.
    w: array of size (N,)
       Contains the weights of the photons that are still propagating, after the absorption event.
    
    
    """ 
    
    indexabs = np.floor(r/dL)
    indexabs = indexabs.astype(int)
    indexabs[indexabs>(Nbins-1)] = Nbins-1
    indexabs[indexabs<0] = 0
    
    indexabs, indices, inverse = np.unique(indexabs,return_index= True, return_inverse=True,return_counts=False,axis=1)
    w_new = np.zeros((np.shape(indexabs)[1],))
    w_new = np.bincount(inverse,w,len(indices))

    A[indexabs[0,:],indexabs[1,:],indexabs[2,:]] += w_new*mu_a/mu_t
    w = w*(1-mu_a/mu_t)
    return A,w


def reflection(r,u,s,L,w,w_refl,n1,n2):
    """reflection(r,u,s,L,w,w_refl,n1,n2):
    Description:
    ------------
    This function reflects photons that are incident at a boundary. A part of the photon weight is transmitted and
    a part of the photon weight is reflected at the boundary, according to the Fresnel equations
    
    Parameters:
    -----------
    r: array of size (3,N)
        array containing the position (x,y,z) of each photon
    u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity component    
    s: array of size (N)
        array containing the distance each photon has travelled
    L   : float
        Length of one side
    w: array of size (N)
        array containing the weight of each photon
    w_refl: float
        sum of the photon weights that are left the volume
    n1: float
        refractive index inside volume
    n2: float
        refractive index inside outside
        
    Results:
    --------
    r: array of size (3,N)
        array containing the position (x,y,z) of each photon after the reflection
    u: array of size (3,N)
        array containing the unit velocity vector after the reflection for all the N photon with the rows containing 
        the x,y,z-velocity component   
    w: array of size (N)
        array containing the weight of each photon    
    w_refl: float
        sum of the photon weights that are left the volume    
    """
    
    N = np.shape(u)
    N = N[1]
    #loop over x-, y- and z-direction
    for i in range(3):
        
        #find photons that are outside the volume
        out_L = r[i,:]>L                #photons outside left boundary
        out_0 = r[i,:]<0                #photons outside right boundary
        out   = out_L+out_0             #photons oustide left and right boundary
        
        #set photons outside volume back to original position
        r[:,out] = r[:,out] - u[:,out]*s[out] 
        
        #Set photons at the boundary of volume
        s1        = np.zeros((N,))
        s1[out_L] = abs(abs(r[i,out_L]-L)/u[i,out_L])
        s1[out_0] = abs(r[i,out_0]/u[i,out_0])
        r[:,out]  = r[:,out] + u[:,out]*s1[out]
        
        #initialize fresnel reflection and transmission
        Rs = np.zeros(np.sum(out))  #s-polarisation reflection
        Rp = np.copy(Rs)            #p-polarisation reflection
        R  = np.copy(Rs)            #total reflection coefficent 
        
        ci = abs(u[i,out])     #cosinus of incident angle
        ci[ci>1] = 1
        si = np.sqrt(1-ci**2)  #sinus of incident angle
       
        tot_ref = n1/n2*si>1        #total internal reflection
        par_ref = tot_ref == False  #reflection and transmission
        
        ci = ci[par_ref]
        si = si[par_ref]
        
        #calculate fresnel coefficients
        Rs[par_ref] = ((n1*ci - n2*np.sqrt(1-(n1/n2*si)**2))/(n1*ci + n2*np.sqrt(1-(n1/n2*si)**2)))**2  #s-polarisation reflection
        Rp[par_ref] = ((-n2*ci + n1*np.sqrt(1-(n1/n2*si)**2))/(n2*ci + n1*np.sqrt(1-(n1/n2*si)**2)))**2 #p-polarisation reflection
        R[par_ref]  = (Rs[par_ref]+Rp[par_ref])/2 
        R[tot_ref]  = 1
                
        #Total change in wheight due to transmitted photons
        w_refl += np.sum((1-R)*w[out])
        
        #update weights, velocity vector
        w[out]   = R*w[out]
        u[i,out] = -u[i,out]
        
        #reflect photons at the boundary
        r[:,out] = r[:,out] + (s[out]-s1[out])*u[:,out]
        
    return r,u,w,w_refl

def terminate(w,u,r,threshold,p_term):
    """terminate(w,threshold,p_term,u,r):
    Description:
    ------------
    This function if a photon is terminated when the weight of the photon drops below a certain threshold value.
    If the weight is below the threshold value then the photon is terminated with a probability of p = 1-p_term, or 
    the weight of the photon is scaled with a factor of 1/p_term with a probability of p = p_term.
    This method of terminating and scaling of the weight ensures that the system conserves energy.
    The photons that are terminated are removed from the input arrays
    
    Parameters:
    -----------
    w: array of size (N)
        array containing the weight of the photons
    u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity 
        component   
    r: array of size (3,N)
        array containing the position (x,y,z) of each photon
    threshold: float
        threshold value if w<threshold than there is a probability that the photon is terminated
    p_term: float
        probability that a photon continues to propagate if w<threshold (i.e 1-p_term is probability that a photon is 
        terminated)
    
    Results:
    --------
    w: array of size (M)
        array containing the weight of the photons with M<=N
    u: array of size (3,M)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity 
        component with M<=N
    r: array of size (3,M)
        array containing the position (x,y,z) of each photon with M<=N
    """
    
    #find weights which are below the treshold value
    thresh = w<threshold
    
    #Deterime which photons will be terminated 
    RND       = np.random.rand(np.sum(thresh),)
    w[thresh] = (RND<p_term) * (1/p_term)*w[thresh]
    live      =  w != 0
    
    #terminate photons
    u = u[:,live]
    r = r[:,live]
    w = w[live]
    return w,u,r


def hop(r,u,mu_t):
    """hop(r,u,mu_t):
    Description:
    ------------
    This function moves photons a random distance from there initial position r along there unit vector u. The each 
    photon travels a random distance s with a probability: p(s) = exp(-mu_t*s)/mu_t
    
    Parameters:
    -----------
    r: array of size (3,N)
        array containing the position (x,y,z) of each photon
    u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity component 
    mu_t: float
        attenuation coefficient
    
    Results:
    --------
    r: array of size (3,N)
        array containing the position (x,y,z) of each photon
    s: array of size (N)
        array containing the distance each photon has travelled    
    """
    
    N = np.shape(r)
    N = N[1]
    RND = np.random.rand(N,)
    s = -np.log(RND)/mu_t
    r = r+u*s

    return r,s


def through_cavity_boundary(r_old,r_new,u,s,w,cavcent,cavr):
    """through_cavity_boundary(r_old,r_new,u,s,w,cavcent,cavr):
    Description:
    ------------
    This function checks which photons cross a cavity boundary and returns the positions, velocities, step sizes, weights etc of these photons and their index in the total array of photons by calculating the analytical solution to the intersection between a line and a sphere.
    
    Parameters:
    -----------
    r_old: array of size (3,N)
        array containing the position (x,y,z) of each photon before the "hop-step"
    r_new: array of size (3,N)
        array containing the position (x,y,z) of each photon after the "hop-step"   
    u: array of size (3,N)
        array containing the unit velocity vector for all the N photon with the rows containing the x,y,z-velocity component    
    s: array of size (N)
        array containing the distance each photon has travelled
    w: array of size (N)
        array containing the weight of each photon
    cavcent: array of size (3)
        array containing the (x,y,z) of a spherical cavity
    cavr: float
        float containing the radius of the spherical cavity
        
    Results:
    --------
    r_intersec: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of the intersection for each photon which path intersects the cavity
    u_intersec: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon which path intersects the cavity
    w_intersec: array of size (number of photons intersecting the cavity)
        array containing the weights of each photon which path intersects the cavity
    dist_to_travel: array of size (number of photons intersecting the cavity)
        array containing the distance each photon has to travel after the intersection with the cavity
    intersec: boolean array of size (N)
        boolean array containing the indices of which photons actually intersect the cavity
    reflec_io: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the inside to the outside of the cavity
    reflec_oi: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the outside to the inside of the cavity
    
    """

    incavity_b = np.sum((r_old-cavcent)**2,axis=0)<cavr**2  #photons before hop in cavity
    incavity_a = np.sum((r_new-cavcent)**2,axis=0)<cavr**2  #photons after hop in cavity

    
    #Only the ones who switch from in to out, or from out to in should be reflected.
    reflec_cav = incavity_b*(np.invert(incavity_a)) + incavity_a*(np.invert(incavity_b)) #photons through cavity
    reflec_io  = incavity_b[reflec_cav]*(np.invert(incavity_a[reflec_cav]))              #photons form inside to outside
    reflec_oi  = incavity_a[reflec_cav]*(np.invert(incavity_b[reflec_cav]))              #photons from outside to insdie
    
    #Find intersections
    #A description of the method used to find the intersection can be found in:
    #http://www.ambrsoft.com/TrigoCalc/Sphere/SpherLineIntersection_.htm
    a = (r_new[0,:]-r_old[0,:])**2 + (r_new[1,:]-r_old[1,:])**2 + (r_new[2,:]-r_old[2,:])**2
    b = 2*((r_new[0,:]-r_old[0,:])*(r_old[0,:]-cavcent[0]) + (r_new[1,:]-r_old[1,:])*(r_old[1,:]-cavcent[1]) + (r_new[2,:]-r_old[2,:])*(r_old[2,:]-cavcent[2]))
    c = np.sum(cavcent**2) + np.sum(r_old**2,axis = 0) - 2*(np.sum(cavcent*r_old,axis = 0)) - cavr**2
    
    D = b**2 - 4*a*c
    intersec = D>0
    intersec = intersec*reflec_cav
    
    t1 = (-b[intersec]+np.sqrt(D[intersec]))/(2*a[intersec]) #parameter of intersection i.e: r = r1+(r2-r1)*t
    t2 = (-b[intersec]-np.sqrt(D[intersec]))/(2*a[intersec]) #parameter of intersection i.e: r = r1+(r2-r1)*t
    
    r1 = r_old[:,intersec] + (r_new[:,intersec]-r_old[:,intersec])*t1              #intersection of trajectory with sphere
    r2 = r_old[:,intersec] + (r_new[:,intersec]-r_old[:,intersec])*t2              #intersection of trajectory with sphere
    
    #calculate nearest intersection coordinate with sphere 
    r1closer = np.sum((r_old[:,intersec] - r1)**2,axis = 0)< np.sum((r_old[:,intersec] - r2)**2, axis = 0) 
    
    #parameters of the photons which go trough the sphere 
    r_intersec = r1closer*r1 + np.invert(r1closer)*r2  #coordinates of the intersection with the sphere
    r_intersec[:,reflec_io] = r1closer[reflec_io]*r2[:,reflec_io] + np.invert(r1closer[reflec_io])*r1[:,reflec_io] 
    u_intersec = u[:,intersec]                         #velocity vectors 
    w_intersec = w[intersec]                           #weight of photons
    dist_to_travel =  np.sqrt(np.sum((r_new[:,intersec] - r_intersec)**2,axis = 0)) #distance photons still have to travel after interaction with sphere
    
    return r_intersec,u_intersec,w_intersec,dist_to_travel,intersec,reflec_io,reflec_oi

def reflect_cavity(r,u,w,intersec,reflec_io,reflec_oi,n1,n2,cavcent):
    """reflect_cavity(r,u,w,dist_to_travel,intersec,reflec_io,reflec_oi,n1,n2,cavcent):
    Description:
    ------------
    This function reflects part of the weight of photons which intersect from the cavity boundary by calculating the Fresnel reflection coefficients (and taking into account total internal reflection).
    
    Parameters:
    -----------
    r: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of the intersection for each photon which path intersects the cavity
    u: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon which path intersects the cavity
    w: array of size (number of photons intersecting the cavity)
        array containing the weights of each photon which path intersects the cavity
    dist_to_travel: array of size (number of photons intersecting the cavity)
        array containing the distance each photon has to travel after the intersection with the cavity
    intersec: boolean array of size (N)
        boolean array containing the indices of which photons actually intersect the cavity
    reflec_io: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the inside to the outside of the cavity
    reflec_oi: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the outside to the inside of the cavity
    n1: float
        float indicating the refraction coefficient of the tissue
    n2: float
        float indicating the refraction coefficient of the cavity
    cavcent: array of size (3)
        array containing the (x,y,z) of a spherical cavity

        
    Results:
    --------
    r: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of each photon which path intersects the cavity
    u: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon which path intersects the cavity
    w: array of size (number of photons intersecting the cavity)
        array containing the weights of each photon which path intersects the cavity
    R: array of size (number of photons intersecting the cavity)
        array containing the fresnel reflection coefficient for each intersecting photon
    ci: array of size (number of photons intersecting the cavity)
        array containing the cosine of the angle between the incoming photons and the cavity surface
    si: array of size (number of photons intersecting the cavity)
        array containing the sine of the angle between the incoming photons and the cavity surface
    normal: array of size (3,number of photons intersecting the cavity)
        array containing vector normal to the cavity surface at the intersection location
    """    
    #Fresnel reflection and transmission
    Rs = np.zeros(np.sum(intersec))   #s-polarisation reflection
    Rp = np.copy(Rs)                  #p-polarisation reflection
    R  = np.copy(Rs)                  #total reflection coefficent 
    
    normal = (r - cavcent)/np.sqrt(np.sum((r - cavcent)**2,axis = 0)) #normal vector sphere
    ci = abs(np.sum(normal*u,axis = 0))                               #cosinus of incident angle
    ci[ci>1] = 1
    si = np.sqrt(1-ci**2)                                             #sinus of incident angle 
    
    #refelction from outside to inside
    tot_ref = n1/n2*si>1        #total internal reflection
    par_ref = tot_ref == False  #reflection and transmission
    oi      = par_ref*reflec_oi 
    
    ci_oi = ci[oi]
    si_oi = si[oi]
    
    #Fresnel coefficents for reflection from outside to inside
    Rs[oi] = ((n1*ci_oi - n2*np.sqrt(1-(n1/n2*si_oi)**2))/(n1*ci_oi + n2*np.sqrt(1-(n1/n2*si_oi)**2)))**2  #s-polarisation
    Rp[oi] = ((-n2*ci_oi + n1*np.sqrt(1-(n1/n2*si_oi)**2))/(n2*ci_oi + n1*np.sqrt(1-(n1/n2*si_oi)**2)))**2 #p-polarisation
    R[oi]  = (Rs[oi]+Rp[oi])/2 
    R[tot_ref*reflec_oi]  = 1
    
    #refelction from inside to outside
    tot_ref = n2/n1*si>1        #total internal reflection
    par_ref = tot_ref == False  #reflection and transmission
    io      = par_ref*reflec_io
    
    ci_io = ci[io]
    si_io = si[io]

    #Fresnel coefficents for reflection from inside to outside
    Rs[io] = ((n2*ci_io - n1*np.sqrt(1-(n2/n1*si_io)**2))/(n2*ci_io + n1*np.sqrt(1-(n2/n1*si_io)**2)))**2  #s-polarisation 
    Rp[io] = ((-n1*ci_io + n2*np.sqrt(1-(n2/n1*si_io)**2))/(n1*ci_io + n2*np.sqrt(1-(n2/n1*si_io)**2)))**2 #p-polarisation 
    R[io]  = (Rs[io]+Rp[io])/2 
    R[tot_ref*reflec_io]  = 1
    
    #calculate new weights of photons and new velocity vectors
    w = R*w
    u = u-2*np.sum(u*normal,axis = 0)*normal 
    
    return u,w,R,ci,si,normal

def transmission_cavity(r,u,w,s,dist_to_travel,ci,si,normal,R,reflec_oi,reflec_io,n1,n2):
    """transmission_cavity(r,u,w,s,dist_to_travel,ci,si,normal,R,reflec_oi,reflec_io,n1,n2):
    Description:
    ------------
    This function transmits and refracts part of the weight of photons which intersect from the cavity boundary by using the calculated Fresnel coefficients, Snell's law and a rotation matrix
    
    Parameters:
    -----------
    r: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of the intersection for each photon which path intersects the cavity
    u: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon which path intersects the cavity
    w: array of size (number of photons intersecting the cavity)
        array containing the weights of each photon which path intersects the cavity
    s: array of size (number of photons intersecting the cavity)
        array containing the step size of each photon which path intersects the cavity
    dist_to_travel: array of size (number of photons intersecting the cavity)
        array containing the distance each photon has to travel after the intersection with the cavity
    ci: array of size (number of photons intersecting the cavity)
        array containing the cosine of the angle between the incoming photons and the cavity surface
    si: array of size (number of photons intersecting the cavity)
        array containing the sine of the angle between the incoming photons and the cavity surface
    normal: array of size (3,number of photons intersecting the cavity)
        array containing vector normal to the cavity surface at the intersection location
    reflec_io: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the inside to the outside of the cavity
    reflec_oi: boolean array of size (N)
        boolean array containing the indices of which photons intersect from the outside to the inside of the cavity
    n1: float
        float indicating the refraction coefficient of the tissue
    n2: float
        float indicating the refraction coefficient of the cavity
        
    Results:
    --------
    utrans: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon which path intersects the cavity and is transmitted and refracted through the surface of the cavity
    rtrans: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of each photon which path intersects the cavity and is transmitted and refracted through the surface of the cavity
    wtrans: array of size (number of photons intersecting the cavity)
        array containing the weights of each photon which path intersects the cavity and is transmitted and refracted through the surface of the cavity
    strans: array of size (number of photons intersecting the cavity)
        array containing the step size of each photon which path intersects the cavity and is transmitted and refracted through the surface of the cavity
    
    """
    
    #find the photons that are partially reflected
    if n1<n2:
        par_ref = reflec_oi+(reflec_io*((n2/n1*si)<1))
    else:
        par_ref = reflec_io+(reflec_oi*((n1/n2*si)<1))
    N = np.sum(par_ref)
    
    #Calculate normal vector of plane of incident for each photon
    axis_r = np.cross(u[:,par_ref].transpose(),normal[:,par_ref].transpose())
    axis_r = axis_r.transpose()
    axis_r = axis_r/np.sqrt(np.sum(axis_r**2,axis = 0))
    axis_r[np.isnan(axis_r)] = 1
    
    
    theta_i = np.arccos(ci[par_ref])      #angle of incident
    theta_u = np.zeros((np.sum(par_ref))) #angle of reflection
    
    theta_u[reflec_oi[par_ref]] = np.arcsin(n1/n2*si[par_ref*reflec_oi])
    theta_u[reflec_io[par_ref]] = np.arcsin(n2/n1*si[par_ref*reflec_io])            #changing angles for which the photons go from in to outside
    
    #rotation angle around axis_r
    theta_rot = theta_u-theta_i
    costheta = np.cos(theta_rot)
    sintheta = np.sin(theta_rot)
               
    #rotation matrix
    Rot_mat = np.array([[costheta + axis_r[0,:]**2*(1-costheta), axis_r[0,:]*axis_r[1,:]*(1-costheta)-axis_r[2,:]*sintheta,                  axis_r[0,:]*axis_r[2,:]*(1-costheta) + axis_r[1,:]*sintheta],
                  [axis_r[0,:]*axis_r[1,:]*(1-costheta)+axis_r[2,:]*sintheta,  costheta + axis_r[1,:]**2*(1-costheta),                 axis_r[1,:]*axis_r[2,:]*(1-costheta) - axis_r[0,:]*sintheta],
                  [axis_r[0,:]*axis_r[2,:]*(1-costheta)-axis_r[1,:]*sintheta,  axis_r[1,:]*axis_r[2,:]*(1-costheta)+axis_r[0,:]*sintheta,  axis_r[2,:]**2*(1-costheta) + costheta]])
               
    #perform rotation
    Rot_mat = np.transpose(Rot_mat,(2,0,1))
    u_rot = u[:,par_ref]
    u_rot = np.transpose(u_rot)
    u_rot = np.reshape(u_rot,((N,3,1)))
    utrans = np.matmul(Rot_mat,u_rot)
    utrans = np.squeeze(np.transpose(utrans,(1,2,0)),axis = 1)
    
    #calculate new weights and distance to travel for the photons 
    wtrans = (1-R[par_ref])*w[par_ref]
    strans = s[par_ref]              
    return utrans,wtrans,strans,par_ref

def cavity_trans_refl(r,r_old,u,w,s,cavcent,cavr,n1,n2,mu_t):
    """cavity_trans_refl(r,r_old,u,w,s,cavcent,cavr,n1,n2,mu_t):
    Description:
    ------------
    This function finds the photons that will traveled trough the cavity. For each photon it determines it position
    after the photon has interacted with the cavity. For each incident photon a part is transmitted, and a part is
    reflected. 
    
    The photons that are transmitted from inside the cavity to the outside of the cavity will travel a random distance 
    according to the mu_t
    
    Parameters:
    -----------
    r: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of each photon after the "hop-step"
    r_old: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of each photon before the "hop-step"
    u: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon 
    w: array of size (N)
        array containing the weight of the photons
    s: array of size (N)
        array containing the distance each photon has travelled    
    cavcent: array of size (3)
        array containing the (x,y,z) of a spherical cavity
    cavr: float
        float containing the radius of the spherical cavity
    n1: float
        float indicating the refraction coefficient of the tissue
    n2: float
        float indicating the refraction coefficient of the cavity
    mu_t: float
        attenuation coefficient
        
    Results:
    --------
    r: array of size (3,number of photons intersecting the cavity)
        array containing the position (x,y,z) of each photon after the interaction
    u: array of size (3,number of photons intersecting the cavity)
        array containing the unit velocity vector of each photon after the interaction
    w: array of size (N)
        array containing the weight of the photons after the interaction 
    out: boolean of size (number of photons intersecting the cavity)  
        boolean containing which photons are inside the cavity (False) and which photons are outside the cavity (True) 
    """

    
    #Compute which photons are in a spherical cavity before the hop
    out = np.sum((r_old-cavcent)**2,axis=0)>cavr**2
    s[out==False] = 2*cavr #ensure that the photons inside the cavity will interact with the boundary  
    
    #find the photons that travel through the cavity and set them at the boundary 
    r_intersec,u_intersec,w_intersec,dist_to_travel,intersec,reflec_io,reflec_oi = through_cavity_boundary(r_old,r,u,s,w,cavcent,cavr)
    
    #reflect the photons that are incident on the boundary of the cavity 
    u[:,intersec],w[intersec],R,ci,si,normal = reflect_cavity(r_intersec,u_intersec,w_intersec,intersec,reflec_io,reflec_oi,n1,n2,cavcent)
    
    #Transmit the photons that are incident on the boundary of the cavity 
    utrans,wtrans,strans,par_ref = transmission_cavity(r_intersec,u_intersec,w_intersec,dist_to_travel,s[intersec],ci,si,normal,R,reflec_oi,reflec_io,n1,n2)
    
    #Look at the photons that have interacted with the boundary of the cavity
    ureflect = u[:,intersec]
    rreflect = np.zeros((3,np.sum(intersec)))
    
    #set the photons that are internal reflect in the cavity just away from the cavity boundary
    rreflect[:,reflec_io] = r_intersec[:,reflec_io] + ureflect[:,reflec_io] * 0.01*cavr
    
    #set the photons that are reflected outside the cavity to their new position
    rreflect[:,reflec_oi] = r_intersec[:,reflec_oi] + ureflect[:,reflec_oi] * dist_to_travel[reflec_oi]
    
    #set the reflected photons back in the original array
    r[:,intersec] = rreflect
    
    #set the transmitted photons to their new position
    rtrans = np.zeros((3,np.sum(par_ref)))
    #photons transmitted from outside cavity to inside cavity
    rtrans[:,reflec_oi[par_ref]] = utrans[:,reflec_oi[par_ref]]*0.01*cavr + r_intersec[:,reflec_oi*par_ref]
    #photons transmitted from inside cavity to outside cavity
    rtrans[:,reflec_io[par_ref]],strans[reflec_io[par_ref]] = hop(r_intersec[:,reflec_io*par_ref],utrans[:,reflec_io[par_ref]],mu_t)

    #add the transmitted photons to the original arrays
    u = np.concatenate((u,utrans),axis = 1)
    r = np.concatenate((r,rtrans),axis = 1)
    w = np.concatenate((w,wtrans),axis = 0)
    s = np.concatenate((s,strans),axis = 0)
    
    #find the photons that are not in the cavity after reflections on cavity boundary
    out = np.sum((r-cavcent)**2,axis=0)>cavr**2  
    return r,u,w,s,out
