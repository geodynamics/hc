SPHERICAL HARMONICS FORMAT

All single layer spherical harmonics are in the fully normalized,
physical convention, e.g. Dahlen & Tromp (1998). 


1) Scalar fields with layers (e.g.: hc_assign_density, which calls:
   sh_read_parameters, sh_init_expansion, sh_read_coefficients)

nr_layers
depth (z, in km)
lmax
a00 b00
a10 b10
a11 b11
a20 b20
a21 b21
....



Unformatted file.


2) Vector (velocity) fields with layers


Unformatted file, velocities are expected to be in cm/yr.



INPUT FILE FORMATS

Description of input file formats, subroutines given in parentheses.


1) Viscosity structure (hc_assign_viscosity)

r_i e_i

Unformatted list of radii (radius of layer/Earth radius) and viscosity
(in Pas) values, reads until end of file. Values determine each layer
viscosity upward until the next entry.

