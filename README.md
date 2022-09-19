# analytic_rupture_radiation


This is a very simple collection of functions and corresponding wrappers allowing to simulate the 3D wavefield of a propagating rupture. If you run into troubles using it, which is quite likely, please do not hesitate to ask the author directly. The code is intended to mimick the observations of ultrafast ultrasound rupture imaging as described in: Aichele, Johannes, soumaya latour, Stefan Catheline, and Philippe Roux. 2022. “Dynamic Full-Field Imaging of Rupture Radiation: Material Contrast Governs Source Mechanism.” Earth and Space Science Open Archive. World. http://www.essoar.org/doi/10.1002/essoar.10512101.1 . The simulation is kinematic and analytic. Therefore the 3D homogeneous elastic Green's functions of a double-couple and a unidirectional singular elastic force are convolved with a source wavelet. To simulate rupture, point sources are superimposed in time and space according to a rupture speed profile. Because ultrasound rupture imaging so far is only implemented in a 2D imaging plane, the ruptrue is a line, where the dimension orthogonal to the rupture direction is omitted/implemented as the surface A of the double couple point source. Please be aware that according to the spatial discretization you choose, you might need a lot of RAM/SWAP (64 Gb and more). The code was last tested in Matlab 2020b. Please feel free to use, modify and cite the code or its corresponding publication according to the gnu license. However, do use it under your own responsibility.

How to test-run code. 
Examples are given in the folder Run files. They correspond to the simulations of the pubication preprint http://www.essoar.org/doi/10.1002/essoar.10512101.1 

General procedure.
Choose the Source point.
Choose if you want to correlate a particle velocity or a displacement source function. 
If you want a single source set Rupture.end and Rupture.start to 0 and the speed constant. 
Hence, Only a single source will be active. 

For a moving Source adapt the Parameters until no error is shown. 

elastography: negative values are upwards polarized
simulation: negative are downwards polarized
--> multiply by -1 to make comparable

90 means right pointing source. 
-1 for leftpointing
+1 for rightpointing

use Vp for small timesteps because otherwise some information from the shear wave might be lost resulting in false amplitudes
