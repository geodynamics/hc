SPHERICAL HARMONICS FORMAT

All single layer spherical harmonics are in the fully normalized,
physical convention, e.g. Dahlen & Tromp (1998). 


All spherical harmonics files have an SH_HEADER

	       lmax layer_i zlabel_i nlayer nrset type

lmax: maximum l of expansion
layer_i:      layer number, from 0....nlayer
zlabel_i:     depth label of this layer, in km (positive, from 0: surface
		     to 2871 for CMB) 
nlayer:	      number of layers
nrset:        number of spherical harmonic coefficient sets    
type:	      0 for Rick, 1 for HealPix etc. (will determine only the
		internal representation, external all is physical
		convention) 

1) Scalar fields with layers (e.g.: hc_assign_density, which calls:
   sh_read_parameters, sh_init_expansion, sh_read_coefficients)


SH_HEADER (see above), with nrset == 1
a00 b00
a10 b10
a11 b11
a20 b20
a21 b21
....



Unformatted file.



2) Surface velocity fields

SH_HEADER, e.g: 127 0 0 1 2 0 (nrset == 2 for poloidal and toroidal fields)
a00p b00p a00t b00t
a10p b10p a10t b10t
...



OTHER INPUT FILE FORMATS

1) Viscosity structure (hc_assign_viscosity)

r_i e_i

Unformatted list of radii (radius of layer/Earth radius) and viscosity
(in Pas) values, reads until end of file. Values determine each layer
viscosity upward until the next entry.



OUTPUT FILES

After a regular run, file sol.bin will have the velocity solution
(cm/yr) in binary format. This file can be extracted using
hc_extract_sh_layer (see below). 

File geoid.ab will have the geoid in meters in a spherical harmonic
expansion.



USING THE OUTPUT

1) Extracting spherical harmonics solutions.

   a) extract SH expansion of radial velocity at layer 2

      hc_extract_sh_layer sol.bin 2 1

   b) extract SH expansion of radial velocity at layer 2 and convert
      to spatial expansion

      hc_extract_sh_layer sol.bin 2 1 0 | sh_syn  

   c) extract SH expansion of poloidal and toroidal velocity at layer 5 and convert
      to spatial expansion of v_theta v_phi

      hc_extract_sh_layer sol.bin 5 2 0 | sh_syn

