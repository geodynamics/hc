INSTALLATION

The code should compile on any basic UNIX/Linux system with the GMT
tools and libraries installed. See the Makefile for comments and what
to modify. 

The Makefile assumes that the following environment variables are
predefined: 

1) ARCH: should be the architecture you are compiling on, e.g. set with

		setenv ARCH `uname -m | gawk '{print(tolower($1))}'`

      for example, i686 for Intel


2) GMTHOME: needs to point to the base directory of the GMT installation,
e.g.
	/usr/local/src/GMTdev/GMT3.4.5/
	

3) NETCDFHOME: for the netcdf libraries, should be e.g.

   /usr/local/src/netcdf-3.5.0/

Also, the Makefile uses the commonly defined compiler variables CC,
F77, CFLAGS, LD, and LDFLAGS.

One common problem is the interface of the ggrd_ subroutines with
GMT. By default, the Makefile is set up for the old syntac of GMT
version 3.4.5. The other alternative is to use 4.1.2 and upward (note
that there were changes between 4.0 and 4.1.2). If you have GMT newer
or equal to 4.1.2, uncomment the two corresponding lines in the
Makefile.

With all things set up, you should be able to type 

make all

to compile the programs.



HC CODE usage

Described in the help page that is displayed for "hc -h" as below. See
also the benchmark scripts in cig_bench.tgz. This tar file will expand
into a cig subdirectory particularly "run_benchmark" there shows how
to run all CIG benchmark tests.


>>>


hc - perform Hager & O'Connell flow computation

This code can compute velocities, tractions, and geoid for simple
density distributions and plate velocities using the semi-analytical
approach of Hager & O'Connell (1981).  This particular implementation
illustrates one possible way to combine the HC solver routines.

Based on code by Brad Hager, Richard O'Connell, and Bernhard Steinberger.

density anomaly options:
-dens   name    use name as a SH density anomaly model (dens.sh.dat)
                This expects density anomalies in % of PREM in DT format.
-ds             additional density anomaly scaling factor (0.2)

Earth model options:
-prem   name    set Earth model to name (/home/becker/progs/src/hc-svn/prem/prem.dat)
-vf     name    set viscosity structure filename to name (visc.dat)
                This file is in non_dim_radius viscosity[Pas] format

boundary condition options:
-fs             perform free slip computation, else no slip or plates (OFF)
-ns             perform no slip computation, else try to read plate velocities (OFF)
-pvel   name    set surface velocity file to name (pvel.sh.dat)
                This file is based on a DT expansion of cm/yr velocity fields.

solution procedure and I/O options:
-ng             do not compute and print the geoid (OFF)
-pptsol         print pol[6] and tor[2] solution vectors (OFF)
-px             print the spatial solution to file (OFF)
-str            compute srr,srt,srp tractions [MPa] instead of velocities [cm/yr] (default: vel)
-v      -vv     -vvv: verbosity levels (0)

<<<

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






DIRECTORY LISTING:

becker@jackie:~/progs/src/hc-svn > ls *.c *.h
ggrd_grdtrack_util.c   hc.h             prem.h          sh_rick.h
ggrd_grdtrack_util.h   hc_init.c        prem_util.c     sh_shana.h
ggrd.h                 hc_input.c       rick_fft_c.c    sh_spherepack.h
ggrd_readgrds.c        hc_matrix.c      rick_sh_c.c     sh_syn.c
ggrd_struc.h           hc_misc.c        sh_ana.c        sh_test.c
ggrd_test.c            hc_output.c      shana_sh.c      simple_test.c
ggrd_velinterpol.c     hc_polsol.c      sh_exp.c        spherepack_sh.c
hc_auto_proto.h        hc_propagator.c  sh.h            test_fft.c
hc_auto_proto.temp.h   hc_solve.c       sh_model.c
hc_extract_sh_layer.c  hc_torsol.c      sh_power.c
hc_filenames.h         main.c           sh_rick_ftrn.h

