README FOR HC

Version 1.0.7 - 2017

Thorsten Becker - twb@ig.utexas.edu


HC is a global mantle circulation solver following Hager & O'Connell
(1981) which can compute velocities, tractions, and geoid for simple
density distributions and plate velocities.

This particular implementation illustrates one possible way to combine
the HC solver routines.

Based on code by Brad Hager, Richard O'Connell, and Bernhard
Steinberger, this version is by Thorsten Becker and Craig O'Neill.

Note also that the Solid Earth Teaching and Research Environment
(SEATREE, http://geosys.usc.edu/projects/seatree/), as preinstalled
with software as a VirtualBox available from the Unified Geodynamics
Earth Science Computing Environment (UGESCE,
http://www-udc.ig.utexas.edu/external/becker/ugesce.html) provides a
convenient user interface for HC.


INSTALLATION

The code should compile on any basic UNIX/Linux system with the GMT
tools (by default, version > 4.2.1, and only version 4.X.X is
supported!) and libraries installed. See the Makefile for comments and
what to modify.

What is described below is the hc Hager & O'Connell (1981) forward
flow and geoid computation. For plate velocity inversions following
the Ricard & Vigny (1989) method, see the hcplates subdirectory and
the README there. The plate inversion code, due to Craig O'Neill isn't
quite working yet.

The Makefile assumes that the following environment variables are
predefined: 

1) GMTHOME: needs to point to the base directory of the GMT installation,
e.g.
	/usr/local/src/GMTdev/GMT4.5.6/

if you installed GMT yourself, or

        /usr/lib/gmt

if you installed GMT with a package manager.

NOTE: ONLY GMT VERSION 4 IS SUPPORTED.

2) NETCDFHOME: for the netcdf libraries.  This could be e.g.

   /usr/local/src/netcdf-3.5.0/

or

   /usr

3) HC_HOME: By default, the object files and binaries will be
   installed in the "objects/" and "bin/" directories in the current
   directory.  If HC_HOME is set (e.g. /usr/local/), then they will
   be installed in $HC_HOME/objects and $HC_HOME/bin.

4) ARCH: If you like, you can put the object files and binaries in
   different directories depending on your architecture.  This may be
   useful if your directories are NFS mounted on different machines.
   One way is to use the output of uname

     setenv ARCH `uname -m | gawk '{print(tolower($1))}'`      # sh/tcsh
     export ARCH=`uname -m | gawk '{print(tolower($1))}'`      # bash

   On a 32 bit Intel machine, this will put the binaries in bin/i686.

Also, the Makefile uses the commonly defined compiler variables CC,
F90, CFLAGS, LD, and LDFLAGS.  So to make static executables, set
LDFLAGS="-static"

By default, the Makefile is set up for the new syntax of GMT version
4.1.2 and higher. The alternative is to use GMT3, which can be done by
defining -DUSE_GMT3 and modifying the two corresponding lines in the
Makefile.include.

With all things set up, you should be able to type 

make all

to compile the programs.



HC CODE usage

Described in the help page that is displayed for "hc -h" as below.
Also see SEATREE (http://geosys.usc.edu/projects/seatree/) for a
graphical user interface, and example plotting scripts as provided
below.

Example input data is provided in subdirectory example_data/

>>>
hc - perform Hager & O'Connell flow computation

This code can compute velocities, tractions, and geoid for simple
density distributions and plate velocities using the semi-analytical
approach of Hager & O'Connell (1981).  This particular implementation
illustrates *one possible way* to combine the HC solver routines,
hc_visc_scan (see below) another. 

Based on code by Brad Hager, Richard O'Connell, and Bernhard
Steinberger.  This version by Thorsten Becker and Craig O'Neill

usage example:

bin/hc -vvv

Compute mantle flow solution using the default input files:
  viscosity profile visc.dat
  density profile   dens.sh.dat
  earth model       prem/prem.dat
and provide lots of output. Default setting is quiet operation.

See README.TXT in the installation directory for example for how to plot output, and
http://geosys.usc.edu/projects/seatree/ for a graphical user interface.
http://www-udc.ig.utexas.edu/external/becker/sdata.html for a VirtualBox install.

density anomaly options:
-dens	name	use name as a SH density anomaly model (dens.sh.dat)
		All density anomalies are in units of 1% of PREM, all SH coefficients
		in Dahlen & Tromp convention.
-dshs		use the short, Becker & Boschi (2002) format for the SH density model (OFF)
-ds		density scaling factor (0.2)
-dnp		do not scale density anomalies with PREM but rather mean density (OFF)
-dsf	file	read depth dependent density scaling from file
		(overrides -ds, OFF), use pdens.py to edit

Earth model options:
-prem	name	set Earth model to name (prem/prem.dat)
-vf	name	viscosity structure filename (visc.dat), use pvisc.py to edit
		This file is in non_dim_radius viscosity[Pas] format
boundary condition options:
-fs		perform free slip computation (ON)
-ns		perform no slip computation (OFF)
-pvel	name	set prescribed surface velocities from file name (OFF)
		The file (e.g. pvel.sh.dat) is based on a DT expansion of cm/yr velocity fields.
-vshs		use the short format (only lmax in header) for the plate velocities (OFF)
-vdir		velocities are given in files name/vel.1.ab to vel.140.ab for different times,
		-140 to -1 Ma before present, where name is from -pvel
-vtime	time	use this particular time step of the plate velocities (-1)

solution procedure and I/O options:
-cbckl	val	will modify CMB boundary condition for all l > val with solver kludge (2147483647)
-ng		do not compute and print the geoid (1)
-ag		compute geoid at all layer depths, as opposed to the surface only
-rg	name	compute correlation of surface geoid with that in file "name",
		this will not print out the geoid file, but only correlations (OFF)
-pptsol		print pol[6] and tor[2] solution vectors (OFF)
-px		print the spatial solution to file (OFF)
-rtrac		compute srr,srt,srp tractions [MPa] instead of velocities [cm/yr] (default: vel)
-htrac		compute stt,stp,spp tractions [MPa] instead of velocities [cm/yr] (default: vel)
-v	-vv	-vvv: verbosity levels (0)



<<<


OTHER BINARIES

hc_visc_scan	Illustrates how to do a parameter space scan for
		viscosities and compute geoid correlations on the
		fly. drive_visc_scan runs a bunch of tests, and
		plot_visc_scan visualizes the output.



KNOWN LIMITATIONS

The propagator matrix approach can be unstable for moderately high
maximum degrees. This behavior can be addressed partially by compiling
HC in quadruple precision, and kinematic boundary can be addressed
with a kludge, see -cbckl. If you compile with this trick and
quadruple precision (the default within SEATREE), you should be good
to up to L=127.


SPHERICAL HARMONICS FORMAT


(A) Regular (long) format, which allows for both scalar and vector
harmonics 
				   
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


(B) Short format

As above, but will only expect a single integer, lmax, as the header
for a spherical harmonic expansion.




OTHER INPUT FILE FORMATS

1) Viscosity structure (hc_assign_viscosity)

r_i e_i

Unformatted list of radii (radius of layer/Earth radius) and viscosity
(in Pas) values, reads until end of file. Values determine each layer
viscosity upward until the next entry. Use the graphical tool
pvisc.py to edit such files.
2) Depth dependent density scaling file 

r_i d_i

Format as for the viscosity file, but d_i are the depth-dependent
scaling factors (this overridings -ds). Use the graphical tool
pdens.py to edit such files. 



OUTPUT FILES

After a regular run, file sol.bin will have the velocity solution
(cm/yr) in binary format. This file can be extracted using
hc_extract_sh_layer (see below). 

File geoid.ab will have the geoid in meters in a spherical harmonic
expansion.



USING THE OUTPUT

1) Extracting spherical harmonics solutions.

   a) extract SH expansion of radial velocity at layer 2

      hc_extract_sh_layer vel.sol.bin 2 1

   b) extract SH expansion of radial velocity at layer 2 and convert
      to spatial expansion

      hc_extract_sh_layer vel.sol.bin 2 1 0 | sh_syn  

   c) extract SH expansion of poloidal and toroidal velocity at layer 5 and convert
      to spatial expansion of v_theta v_phi

      hc_extract_sh_layer vel.sol.bin 5 2 0 | sh_syn



      Also see the script calc_vel_and_plot for some suggestions on
      how to convert the output.

   d) convert velocity output and density anomalies into a binary VTK file
      for paraview

       hc_extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > vel.vtk


       Load a view into paraview (pvsm file is included here)

       paraview --state=hc_world.pvsm

       to get the view in hc_world.png.


USING THE SPHERICAL HARMONICS TOOLS AS STANDALONE


1) Convert scalar values from a GMT/Netcdf grd file into a spherical
   harmonics expansion of degree 31


   a) obtain scalar data at the Gaussian intergration points

      sh_ana -127 | grdtrack -Lx -G$datadir/etopo5/etopo5.1.grd > etopo5.dat


   b) expand

       cat etopo5.dat | sh_ana 127 > etopo5.ab

   c) Alternative: use the grd file directly

   sh_ana 127 $datadir/etopo5/etopo5.1.grd

2) convert spherical harmonics to spatial expansion

      cat etopo5.ab | sh_syn > etopo5.127.dat


Note that sh_syn and sh_ana are only example implementations of the
subroutines, there's very limited actual functionality. For a more
useful spherical harmonics package, see shansyn at 

http://www-udc.ig.utexas.edu/external/becker/sdata.html


That being said, also note helper programs sh_corr and sh_power.


SEATREE

HC is a module of the Solid Earth Teaching and Research Environment
(SEATREE) which provides a graphical user interface to flow
computations and plotting. 

https://geosys.usc.edu/projects/seatree/wiki/ 
     
UGESCE

The Unified Geodynamics Earth Science Computing Environment (UGESCE,
http://www-udc.ig.utexas.edu/external/becker/ugesce.html) provides a VirtualBox
Linux install that includes SEATREE, HC, and a range of other Earth
Science data and software, all in one (big) package, ready to go.




				   

