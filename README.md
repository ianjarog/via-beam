Calculates beam response of the VLA as 
a function of radial distance.

Parameters
----------
inpmap: string
    name of 2D fits file from which coordinates are taken
rad: array
    input radius in arcminute, i.e, the maximum radius  

Methods:
--------
get_vla_pbc()
    return 2D fits file containing the primary beam response values. 

Examples
--------
response = Getpbresponse("fits_file.fits", rad=30, output="vla_pbcmap.fits") #30 arcminute
vla_beam = response.get_vla_pbc()
