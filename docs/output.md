Output
======

You will find the following files after run the tracking code.

vormax\_loc\_\*\*\*\*.txt

vortrack\_\*\*\*\*\_\*\*\*\*.txt

vor\_out\_\*\*\*\*.dat

In vormax\_loc\_(kt).txt, the locations of the vortices at time step =
kt are written in the following format.

|x or longitude  |y or latitude  |Vorticity (10^-3^ s^-1^)  |Size of vortex (km^2^)  |Type
|----------------|---------------|--------------------------|------------------------|------

(Type = 0: mesocyclone, Type = 1: a part of cold front, Type = 2:
synoptic-scale low, Type = 3: synoptic-scale low with cold front)

If several vortices merged into the i-th vortex, the tracks of such
vortices are written in vortrack\_(i)\_0001.txt,
vortrack\_(i)\_0002.txt, vortrack\_(i)\_0003.txt,...

The track of a vortex are written in the following format.

|x or longitude  |y or latitude  |Vorticity (10^-3^ s^-1^)  | Date and time |Size of vortex (km^2^)  |Type
|----------------|---------------|--------------------------|---------------|------------------------|------

Vor\_out\_(kt).dat contains original vorticity field, smoothed vorticity
field, vortex area, x and y component of steering wind, and SLP in grads
format.
