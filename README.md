# Sensor Array Response Simulator

#### A library that simulates array responses for sensors with arbitrary geometry and directional characteristics.

---
>   Archontis Politis, 2014 
>
>   Department of Signal Processing and Acoustics, Aalto University, Finland 
>   archontis.politis@aalto.fi 
---

This is a collection of MATLAB routines for simulation of array responses of

a) directional sensors and,

b) sensors mounted on, or at a distance from, a rigid spherical/cylindrical 
scatterer. 

The computation of their frequency and impulse responses is 
based on the theoretical expansion of a scalar incident plane wave field to 
a series of wavenumber-dependent Bessel-family functions and 
direction-dependent Fourier or Legendre functions.

A function for arbitrary open arrays of directional microphones is included 
not based on the expansion but directly on the steering vector formula of 
inter-sensor delays and sensor gains for directional patterns 
(e.g. arrays of cardioid microphones).

Most of the functionality of the library is displayed at [http://research.spa.aalto.fi/projects/arraysim-lib/arraysim.html], 
or in the included script TEST_ARRAY_SIMULATOR.m

For more information on the expansions, you can have a look on

    Earl G. Williams, "Fourier Acoustics: Sound Radiation and Nearfield 
    Acoustical Holography", Academic Press, 1999

    Heinz Teutsch, "Modal Array Signal Processing: Principles and 
    Applications of Acoustic Wavefield Decomposition", Springer, 2007

and for example on
[http://en.wikipedia.org/wiki/Plane_wave_expansion] and
[http://en.wikipedia.org/wiki/Jacobi-Anger_expansion]

For any questions, comments, corrections, or general feedback, please
contact archontis.politis@aalto.fi

---

For more details on using functions, check their help output in Matlab.

### List of MATLAB files:

* sph_besselj.m     :   Various spherical Bessel-family functions and their
* sph_bessely.m     :   derivatives, for the computation of the radial 
* sph_function.m    :   terms of the expansions (modal weights)
* sph_hankel1.m     
* sph_hankel2.m     
* dbesselj.m        
* dbessely.m        
* dhankel1.m
* dhankel2.m
* dsph_besselj.m
* dsph_bessely.m
* dsph_function.m
* dsph_hankel1.m
* dsph_hankel2.m

* simulateCylArray.m : Simulate cylindrical open arrays or arrays of sensors mounted on a rigid cylinder
* simulateSphArray.m : Simulate spherical open arrays or arrays of sensors mounted on a rigid sphere
* cylModalCoeffs.m : Compute the radial (modal) weights for a cylindrical scatterer
* sphModalCoeffs.m : Compute the radial (modal) weights for a spherical scatterer

* sphericalScatterer.m : Compute the response for a measurement point at an arbitrary distance from a rigid spherical scatterer
* cylindricalScatterer.m : Compute the response for a measurement point at an arbitrary distance from a rigid cylindrical scatterer

* getArrayResponse.m' :  Simulate open arrays of arbitrary directional microphones with axisymmetric responses
