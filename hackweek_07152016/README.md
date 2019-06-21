#Galaxy-galaxy Strong Lensing in PhoSim

### Introduction
In the Twinkles simulation, multiply lensed quasars and supernovas are presented well using the lensing system from OM10 and PhoSim. But the images of lenses and the lensed images of host galaxies of quasars or supernovas are needed. We developed an approach to produce images of lenses and lensed host galaxies, and combine them together using PhoSim according to a given lens from OM10.

### Images of Lenses.  
- Randomly pick one lensing system from OM10 (see Table 1).
- Find the closest velocity dispersion from a galaxy catalog (Table 1 in [paper 1](http://arxiv.org/abs/1501.04977v2)) to obtain Effective radius and Apparent magnitude.
- Assume that the Orientation, Axis ratio, and Redshift of the images are the same with the lens.
- Compose a python script to create an instance catalog based on these parameters.
- Produce the image of the lens in 6 bands by feeding the catalog into PhoSim.

     |OM10                      |  Parameters of the image
:-------:|:-------------------------:|:-------------------------:
ID    | 630751| Same
Orientation  (degree) | 156.5082 | Same
Axis Ratio (a/b) |1.117 | Same
Velocity Dispersion | 287.7398 (km/s) |Same
Redshift of the Lens| 0.228| Same
Effective radius|         -          |1.526 "
Apparent magnitude|     -        | 17.654
Index of Sersic|          -          | 4  
Redshfit of the source| zs = 2.54 | -


### Lensed Images
- Use the same lensing system picked from OM10.
- Randomly generate an instance catalog of a source galaxy.
- Run a ray-tracing simulation to produce lensed arc (fits file).
- Calculate the total magnitude of lensed images according to the input catalog.
- Input the fits file, apparent magnitude, and a reasonable SED into PhoSim.


### Galaxy-galaxy Strong lensing in PhoSim

![ugr](https://cloud.githubusercontent.com/assets/1016652/17046689/b4c857da-4f9b-11e6-910c-ed2ea0552832.jpg)
![gri](https://cloud.githubusercontent.com/assets/1016652/17046690/b6f9755c-4f9b-11e6-8bb2-b4f6946ddb3e.jpg)

### Future work
- Images of Lenses
    - Try Fundamental planes with larger samples of observed galaxies, and more realistic index of Sersic profile. 
    - How to pick a reasonable SED according to the parameters we have so far? (Thomas Collett's Code: [LensPoP](https://github.com/tcollett/LensPop)) 
- Images of Sources 
     - Create the instance catalog of source galaxies (host galaxies) according to the lensed quasars or supernovas.
     - Angular position, ellipticity, orientation, Effective radius, apparent magnitude, SED....
- Combine lens galaxies, lensed galaxies, and lensed quasars (supernovas) together for Twinkles.
- Think about integrating ray-tracing routines into PhoSim to parameterize the strongly lensed arcs.

### References 
- http://arxiv.org/abs/1501.04977v2
- https://github.com/drphilmarshall/OM10
