# Monte Carlo simulation of radiative dose deposition 

The fourth and final project of the computational physics course is about radiative dose deposition. Monte Carlo techniques are used to simulate the path of many 'photons', dropping dose and scattering as they travel.


## Getting Started
In order to run the simulation, copy all the .py files in the repository to your working directory. Change the input constants in constants.py, save the file, then start the simulation from the simulation module, and use the results to do the processing, as follows:

```
import importlib
import simulation
import data_processing

importlib.reload(simulation)
results = simulation.start()
importlib.reload(data_processing)
data_processing.start(*results)
```
Or, simply run launcher.py which does all of the above steps.

### Options for constants

The following inputs can be adjusted in constants.py

Simulation constants
* **Nbins**         The number of bins in each direction. This constant decides what your resolution will be
* **L**             The size of the volume in units of length in which the photon travel will take place.
* **treshold**      The fraction of photon remaining at which it will undergo the "Terminate?" routine. 
* **p_term**        The probability of terminating a photon reaching a weight below threshold value

Material constants
* **mu_a**          Absorption coefficient
* **mu_s**          Scattering coefficient
* **g**             Anisotropy of scattering
* **n1**            Refractive index inside volume
* **n2**            Refractive index outside volume and inside cavity

Beam properties
* **N**             Number of photons to trace
* **a**             Radius of the circular beam
* **mode**          Type of beam. Options are "single" for a single beam pointing in the z-direction, or "multiple"
* **nbps**          Number of beams per side in case "multiple" was chosen for mode

Cavity properties
* **cavity**       Enables the cavity when set to 'True'
* **cavcent**      Coordinates of the center of the cavity
* **cavr**         Radius of the cavity

Tissue properties
* **t_c**          Coordinates of the center of the tumor
* **t_r**          Radius of the tumor
* **spine_c**      Coordinates of the axis of the cylindrical spine
* **spine_r**      Radius of the spine

### Obtained results
The output of the simulation contains:

* **A**            Matrix of Nbins x Nbins x Nbins, containing the energy absorped in each bin, in units of photon weights
* **dose_t**       Total photons weights absorbed in tumor
* **dose_spine**   Total photon weights absorbed in spine


## Authors
* Martijn Eppenga
* Richard Faasse
* Cyrus Tirband 
