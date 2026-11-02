# HC: a global mantle circulation solver

Version 1.1.3, November 2026

Thorsten Becker, twb@ig.utexas.edu, https://twb.ig.utexas.edu/

HC is a global mantle circulation solver following Hager & O'Connell
(1981) which can compute velocities, tractions, and the geoid for
arbitrary density distributions and purely radial viscosity
variations.

This particular implementation illustrates one possible way to combine
the HC solver routines. There is a forward tool example, as included
in the graphical user interface SEATREE, and a parameter space
exploration tool to scan viscosity distributions in terms of their
predictions of the geoid.

This code is based on Fortran routines by Brad Hager, Richard
O'Connell, and most recently Bernhard Steinberger, who made several
modifications. This version is by Thorsten Becker with contributions
to the plate velocity inversion tool by Craig O'Neill (the latter
remains incomplete). The multiple shooting poloidal solver
(`src/hc_polsol_ms.c`), the banded LU routines, the regression test
suite (`test/`) and a number of fixes were written by Claude
(Anthropic) in collaboration with Thorsten Becker in November 2026.

## Contents

- [Directory layout](#directory-layout)
- [How to cite](#how-to-cite)
- [Installation](#installation)
- [Running hc](#running-hc)
- [Other binaries and tools](#other-binaries-and-tools)
- [Solvers and known limitations](#solvers-and-known-limitations)
- [Tests](#tests)
- [File formats](#file-formats)
- [Using the output](#using-the-output)
- [Spherical harmonics tools](#spherical-harmonics-tools)
- [SEATREE and UGESCE](#seatree-and-ugesce)
- [References](#references)
- [Copyright](#copyright)

## Directory layout

| Directory | Content |
|---|---|
| `src/` | C sources and headers for `hc` and the related tools |
| `tools/` | helper scripts (plotting, benchmarks, drivers), see below |
| `test/` | regression tests, run with `make test` |
| `example_data/` | example density, viscosity, plate velocity and reference geoid files |
| `prem/` | PREM model files |
| `tomography/` | a collection of tomographic models in the hc short format |
| `viscosity/` | viscosity profiles |
| `dscale/` | depth dependent density scaling files |
| `geoid/` | observed geoid expansions (EGM2008, ITG) |
| `hcplates/` | plate velocity inversion (incomplete, separate Makefile and README) |
| `anisotropic_viscosity/` | anisotropic viscosity experiments |
| `docker/` | Dockerfile for a container with a built version |
| `bin/`, `objects/` | created by `make` |

## How to cite

If you use HC, please cite the software and the papers that describe
the method, for example as

> Hager, B.H., O'Connell, R.J., Steinberger, B., and Becker, T.W. (2022),
> HC v1.0.15 [software], Computational Infrastructure for Geodynamics,
> https://geodynamics.org/resources/hc

or, for the GitHub version,

> Becker, T.W., O'Neill, C., and Steinberger, B. (2009), HC, a global
> mantle circulation solver, available online at
> https://github.com/geodynamics/hc

together with Hager & O'Connell (1981), Hager & Clayton (1989), and
Steinberger (2000), see [References](#references).

## Installation

The code should compile on any basic UNIX/Linux system. The core
programs need nothing but a C compiler and `make`; the spherical
harmonics tools `sh_syn` and `sh_ana` and the GMT based tools
additionally need the GMT 4 libraries (only GMT 4.x is supported).

### Without GMT

```
make all_no_gmt ADD_FLAGS=-DNO_GMT
```

or switch the include line in the `Makefile` to
`Makefile.include_no_gmt`. This builds `hc`, `hc_visc_scan`,
`hc_invert_dtopo`, `hc_extract_sh_layer`, `hc_extract_spatial`,
`rotvec2vel`, `print_gauss_lat`, `sh_corr`, `sh_power`,
`sh_extract_layer` and `gaussp` into `bin/`.

### With GMT 4

The Makefile assumes that the following environment variables are
set:

1. `GMT4HOME`: the base directory of the GMT 4 installation,
   e.g. `/usr/local/src/GMTdev/GMT4.5.18/` if you installed GMT
   yourself, or `/usr/lib/gmt` if you installed GMT with a package
   manager.
2. `NETCDFHOME`: the base directory of the netcdf libraries,
   e.g. `/usr/local/src/netcdf-3.5.0/` or `/usr`.
3. `HC_HOME` (optional): by default, object files and binaries go to
   `objects/` and `bin/` in the current directory. If `HC_HOME` is
   set (e.g. `/usr/local/`), they go to `$HC_HOME/objects` and
   `$HC_HOME/bin`.
4. `ARCH` (optional): to put object files and binaries into
   architecture dependent subdirectories, e.g. for NFS mounted
   directories shared between machines. One way is to use the output
   of `uname`:
   ```
   export ARCH=`uname -m | gawk '{print(tolower($1))}'`
   ```

The Makefile uses the commonly defined compiler variables `CC`,
`F90`, `CFLAGS`, `LD`, and `LDFLAGS`. To make static executables,
set `LDFLAGS="-static"`. With all things set up,

```
make all
```

compiles all programs. See the `Makefile` and `Makefile.include`
for comments and what to modify. The GMT and netcdf libraries are
only linked into the GMT based tools.

What is described below is the `hc` Hager & O'Connell (1981) forward
flow and geoid computation. For plate velocity inversions following
the Ricard & Vigny (1989) method, see the `hcplates/` subdirectory
and the README there. The plate inversion code is not quite working
yet.

## Running hc

Usage is described in the help page that is displayed for `hc -h`,
reproduced below. The default input file names assume that `hc` is
run from the main hc directory, e.g.

```
bin/hc -vvv
```

computes a flow solution for the example density model and viscosity
structure in `example_data/` with lots of output. Also see SEATREE
(https://github.com/thwbecker/seatree/) for a graphical user
interface, and the example scripts in `tools/`.

```
hc - perform Hager & O'Connell flow computation

This code can compute velocities, tractions, and geoid for simple density distributions
and plate velocities using the semi-analytical approach of Hager & O'Connell (1981).
This particular implementation illustrates one possible way to combine the HC solver routines.
Based on code by Brad Hager, Richard O'Connell, and Bernhard Steinberger.
This version by Thorsten Becker, with contributions by Craig O'Neill.
Copyright (c) 2004-2026 Thorsten Becker, GPL, see README.md
usage example:

hc -vvv

Compute mantle flow solution using the default input files:
  viscosity profile example_data/visc.dat
  density profile   example_data/dens.sh.dat
  earth model       prem/prem.dat
(i.e. hc is expected to be run from the main hc directory).
and provide lots of output. Default setting is quiet operation.

See README.md in the installation directory for example for how to plot output, and
https://github.com/thwbecker/seatree/ for a graphical user interface.
https://twb.ig.utexas.edu/ugesce.html for a VirtualBox install.

density anomaly options:
-dens	name	use name as a SH density anomaly model (example_data/dens.sh.dat)
		All density anomalies are in units of 1% of PREM, all SH coefficients
		in Dahlen & Tromp convention.
-dshs		use the short, Becker & Boschi (2002) format for the SH density model (OFF)
-ds	val	density scaling factor (0.2)
-dnp		do not scale density anomalies with PREM but rather mean density (OFF)
-dsf	file	read depth dependent density scaling from file
		(overrides -ds, OFF), use pdens.py to edit

Earth model options:
-prem	name	set Earth model to name (prem/prem.dat)
-vf	name	viscosity structure filename (example_data/visc.dat), use pvisc.py to edit
		This file is in non_dim_radius viscosity[Pas] format

boundary condition options:
-fs		perform free slip computation (ON)
-ns		perform no slip computation (OFF)
-pvel	name	set prescribed surface velocities from file name (OFF)
		The file (e.g. example_data/pvel.sh.dat) is based on a DT expansion of cm/yr velocity fields.
-vshs		use the short format (only lmax in header) for the plate velocities (OFF)
-rnr		remove any net rotation component from the plate velocities (OFF)
-vdir		velocities are given in files name/vel.1.ab to vel.140.ab for different times,
		-140 to -1 Ma before present, where name is from -pvel
-vtime	time	use this particular time step of the plate velocities (-1)

solution procedure and I/O options:
-solver	val	poloidal solver: 0: original propagator matrix (inaccurate for l > ~50), 1: multiple shooting (stable for high l) (1)
-msg	val	multiple shooting: max growth factor (r2/r1)^(2l+3) per interval (10000)
-cbckl	val	will modify CMB boundary condition for all l > val with solver kludge (2147483647)
-ng		do not compute and print the geoid (1)
-ag		compute geoid at all layer depths, as opposed to the surface only
-pptsol		print pol[6] and tor[2] solution vectors (OFF)
-px		print the spatial solution to file (OFF)
-ascii		write the spherical harmonic solution in ASCII instead of binary (OFF)
-rtrac		compute srr,srt,srp tractions [MPa] instead of velocities [cm/yr] (default: vel)
-htrac		compute stt,stp,spp tractions [MPa] instead of velocities [cm/yr] (default: vel)
-pk		print kernels only, up to default L = 31
-v	-vv	-vvv: verbosity levels (0)
```

## Other binaries and tools

| Binary | Purpose |
|---|---|
| `hc` | forward flow, traction and geoid computation |
| `hc_visc_scan` | parameter space scan over layered viscosities, computes geoid correlations on the fly; `tools/drive_visc_scan` runs a set of scans, `tools/plot_visc_scan` plots them |
| `hc_invert_dtopo` | dynamic topography inversion (experimental) |
| `hc_extract_sh_layer` | extract spherical harmonic expansions of one layer from the binary solution files |
| `hc_extract_spatial` | convert the solution to VTK files for paraview |
| `sh_syn`, `sh_ana` | spherical harmonics synthesis and analysis (GMT build only) |
| `sh_corr`, `sh_power` | correlation and power of spherical harmonic expansions |
| `rotvec2vel` | plate velocities from rotation vectors (NNR-MORVEL56) |

The scripts in `tools/` find the hc directory from their own location,
so they can be run from the main directory as `tools/script`, as long
as the code has been compiled (binaries in `bin/`):

| Script | Purpose |
|---|---|
| `calc_vel_and_plot` | compute velocities or tractions and plot them with GMT |
| `calc_dyn_topo` | compute radial tractions and dynamic topography |
| `compute_and_extract_tractions` | tractions for free slip and plate boundary conditions |
| `plot_geoid` | predicted vs. observed geoid map |
| `plot_spectral_model`, `plot_spatial_model` | plot layers of a solution |
| `hc2gmt`, `calc_rms` | convert spatial output (`-px`) to GMT grids and compute RMS values |
| `drive_visc_scan`, `plot_visc_scan` | viscosity scans and plots |
| `pvisc.py`, `pdens.py` | graphical editors for viscosity and density scaling files (python 2, GTK) |
| `test_sh_routines`, `test_codes`, `test_grdgrd2correlation`, `compare_solution_methods`, `run_bench` | legacy tests and benchmarks, partly relying on tools outside of hc |
| `make_binary_package` | build a binary distribution directory |

Several of the plotting scripts need GMT, and some rely on further
tools (e.g. `shsyn` from the shansyn package) that are noted in the
script headers.

## Solvers and known limitations

The original, propagator matrix solution of the poloidal problem
(`src/hc_polsol.c`, `-solver 0`) becomes inaccurate for moderately
high maximum degrees (roughly L > 50 in double precision, depending
on the viscosity structure; at L = 127, degrees above about 80 are
wrong by order unity). The CMB boundary condition can be modified for
high degrees with a kludge following Steinberger & Torsvik (2012),
see `-cbckl`, which only applies to this solver.

Since November 2026, the default is a multiple shooting solver,
`-solver 1` (`src/hc_polsol_ms.c`, written by Claude/Anthropic with
T. Becker). It uses the same propagators, sources, and boundary
conditions, but solves for the state vectors at all levels at once in
a banded linear system, which stays accurate for arbitrary L in
double precision. It agrees with the propagator solution to about
1e-10 where the latter is accurate, is insensitive to the layering,
and has been tested up to L = 255. The levels are subdivided such
that the propagator growth per interval, (r2/r1)^(2L+3), stays below
the `-msg` value (1e4 by default; results do not depend on this
setting in practice). The cost is about one second for L = 255 with
26 layers. The header of `src/hc_polsol_ms.c` explains the method and
the reasons for the instability of the propagator approach in more
detail. Earlier versions offered a long double ("quad precision")
compile option; this was removed since it did not compile with the
generated prototypes and only addressed part of the problem.

The plate velocity inversion (`hcplates/`) and the dynamic topography
inversion (`hc_invert_dtopo`) are incomplete.

## Tests

| Target | What it does |
|---|---|
| `make test` | runs the regression tests in `test/` (using `example_data/`) for both the propagator matrix (`-solver 0`) and the multiple shooting (`-solver 1`) solver against reference output computed with the propagator solver at lmax = 31, and checks the stability of the multiple shooting solver at lmax = 127 |
| `make test_asan` | rebuilds with address/undefined behavior sanitizers and runs the tests with those binaries |

`test/propagator_stability.sh` is a diagnostic that measures the
round off amplification of a solver with increasing lmax.

## File formats

### Spherical harmonics

All single layer spherical harmonics are in the fully normalized,
physical convention, e.g. Dahlen & Tromp (1998).

**(A) Regular (long) format**, which allows for both scalar and vector
harmonics. All spherical harmonics files have a header

```
lmax layer_i zlabel_i nlayer nrset type
```

| Field | Meaning |
|---|---|
| `lmax` | maximum l of expansion |
| `layer_i` | layer number, from 0 ... nlayer |
| `zlabel_i` | depth label of this layer, in km (positive, from 0 at the surface to 2871 at the CMB) |
| `nlayer` | number of layers |
| `nrset` | number of spherical harmonic coefficient sets |
| `type` | 0 for Rick, 1 for HealPix etc. (determines only the internal representation, externally all is physical convention) |

1. Scalar fields with layers (e.g. `hc_assign_density`, which calls
   `sh_read_parameters`, `sh_init_expansion`, `sh_read_coefficients`):
   header with `nrset == 1`, followed by the coefficients, unformatted
   ```
   a00 b00
   a10 b10
   a11 b11
   a20 b20
   a21 b21
   ...
   ```
2. Surface velocity fields: header, e.g. `127 0 0 1 2 0` (`nrset == 2`
   for poloidal and toroidal fields), followed by
   ```
   a00p b00p a00t b00t
   a10p b10p a10t b10t
   ...
   ```

**(B) Short format**, as above, but with only a single integer,
`lmax`, as the header of each spherical harmonic expansion (the
tomographic models in `tomography/` use this format, select it with
`-dshs`).

### Other input files

1. Viscosity structure (`hc_assign_viscosity`, `-vf`): an unformatted
   list of
   ```
   r_i e_i
   ```
   with radius (radius of layer / Earth radius) and viscosity (in Pas),
   read until end of file. Each value determines the viscosity upward
   until the next entry. Use the graphical tool `tools/pvisc.py` to
   edit such files.
2. Depth dependent density scaling (`-dsf`): same format, but `d_i`
   are the depth dependent scaling factors (this overrides `-ds`).
   Use `tools/pdens.py` to edit such files.

### Output files

After a regular run, `vel.sol.bin` holds the velocity solution (cm/yr)
in binary format (`rtrac.sol.bin` for tractions with `-rtrac`,
`vel.sol.dat` in ASCII with `-ascii`). These files can be read with
`hc_extract_sh_layer`. `geoid.ab` holds the geoid in meters as a
spherical harmonic expansion, `vdepth.dat` the depths of the solution
layers, and `dscaled.sol.bin` the scaled density anomalies.

## Using the output

1. Extract the spherical harmonic expansion of the radial velocity at
   layer 2:
   ```
   hc_extract_sh_layer vel.sol.bin 2 1
   ```
2. The same, converted to a spatial expansion:
   ```
   hc_extract_sh_layer vel.sol.bin 2 1 0 | sh_syn
   ```
3. Poloidal and toroidal velocity at layer 5, converted to a spatial
   expansion of v_theta, v_phi:
   ```
   hc_extract_sh_layer vel.sol.bin 5 2 0 | sh_syn
   ```
   Also see `tools/calc_vel_and_plot` for suggestions on how to
   convert and plot the output.
4. Convert velocity output and density anomalies into a binary VTK
   file for paraview, and load the view in `tools/hc_world.pvsm` to
   get `tools/hc_world.png`:
   ```
   hc_extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > vel.vtk
   paraview --state=tools/hc_world.pvsm
   ```

## Spherical harmonics tools

1. Convert scalar values from a GMT/netcdf grd file into a spherical
   harmonics expansion of degree 127: obtain the data at the Gaussian
   integration points, then expand
   ```
   sh_ana -127 | grdtrack -Lx -G$datadir/etopo5/etopo5.1.grd > etopo5.dat
   cat etopo5.dat | sh_ana 127 > etopo5.ab
   ```
   or use the grd file directly
   ```
   sh_ana 127 $datadir/etopo5/etopo5.1.grd
   ```
2. Convert spherical harmonics to a spatial expansion:
   ```
   cat etopo5.ab | sh_syn > etopo5.127.dat
   ```

Note that `sh_syn` and `sh_ana` are only example implementations of
the subroutines, with very limited functionality. For a more useful
spherical harmonics package, see shansyn at
https://twb.ig.utexas.edu/sdata.html. Also note the helper programs
`sh_corr` and `sh_power`.

## SEATREE and UGESCE

HC is a module of the Solid Earth Teaching and Research Environment
(SEATREE, https://github.com/thwbecker/seatree), which provides a
graphical user interface to flow computations and plotting, for
selecting seismic tomography models, editing viscosity structures,
etc.

The Unified Geodynamics Earth Science Computing Environment (UGESCE,
https://twb.ig.utexas.edu/ugesce.html) provides a VirtualBox Linux
install that includes SEATREE, HC, and a range of other Earth Science
data and software, all in one (big) package, ready to go.

## References

- Becker, T.W., and Boschi, L. (2002), A comparison of tomographic and
  geodynamic mantle models, Geochem. Geophys. Geosyst., 3(1), 1003,
  doi:10.1029/2001GC000168
- Dahlen, F.A., and Tromp, J. (1998), Theoretical Global Seismology,
  Princeton University Press
- Dziewonski, A.M., and Anderson, D.L. (1981), Preliminary reference
  Earth model, Phys. Earth Planet. Inter., 25, 297-356
- Hager, B.H. (1984), Subducted slabs and the geoid: constraints on
  mantle rheology and flow, J. Geophys. Res., 89, 6003-6015
- Hager, B.H., and Clayton, R.W. (1989), Constraints on the structure
  of mantle convection using seismic observations, flow models, and
  the geoid, in: Mantle Convection, edited by W.R. Peltier, Gordon and
  Breach, New York, 657-763
- Hager, B.H., and O'Connell, R.J. (1981), A simple global model of
  plate dynamics and mantle convection, J. Geophys. Res., 86,
  4843-4867
- Panasyuk, S.V., Hager, B.H., and Forte, A.M. (1996), Understanding the
  effects of mantle compressibility on geoid kernels, Geophys. J. Int.,
  124, 121-133
- Ricard, Y., and Vigny, C. (1989), Mantle dynamics with induced plate
  tectonics, J. Geophys. Res., 94, 17543-17559
- Richards, M.A., and Hager, B.H. (1984), Geoid anomalies in a dynamic
  Earth, J. Geophys. Res., 89, 5987-6002
- Steinberger, B. (2000), Slabs in the lower mantle - results of dynamic
  modelling compared with tomographic images and the geoid,
  Phys. Earth Planet. Inter., 118, 241-257
- Steinberger, B., and Torsvik, T.H. (2012), A geodynamic model of plumes
  from the margins of Large Low Shear Velocity Provinces,
  Geochem. Geophys. Geosyst., 13, Q01W09, doi:10.1029/2011GC003808

Other tools in this distribution:

- Argus, D.F., Gordon, R.G., and DeMets, C. (2011), Geologically current
  motion of 56 plates relative to the no-net-rotation reference frame,
  Geochem. Geophys. Geosyst., 12, Q11001, doi:10.1029/2011GC003751
  (NNR-MORVEL56 rotation vectors in `rotvec2vel`)
- Milner, K., Becker, T.W., Boschi, L., Sain, J., Schorlemmer, D., and
  Waterhouse, H. (2009), The Solid Earth Research and Teaching
  Environment: a new software framework to share research tools in
  the classroom and across disciplines, Eos Trans. AGU, 90(12), 104
  (SEATREE)

## Changes

- 11/2026: multiple shooting poloidal solver (`src/hc_polsol_ms.c`,
  Claude/Anthropic with T. Becker) added and made the default
  (`-solver 1`); the original propagator matrix solver remains
  available as `-solver 0`. Banded LU routines in `src/hc_matrix.c`.
  Regression test suite in `test/` (`make test`, `make test_asan`).
  Sources moved to `src/`, scripts to `tools/` (they locate the hc
  directory themselves and can be run as `tools/script` from the main
  directory), data files to `example_data/`, `dscale/`, `geoid/`,
  `viscosity/`; default input file names now point to
  `example_data/`. README converted to Markdown. Fixes: heap overflow
  in the four layer viscosity mode (`hc_visc_scan`), out of bounds
  read in the `dvisc` assignment, `nprops_max` too small for few
  density layers, correlation lmax typo, 3 x 3 solver rewritten (no
  longer Numerical Recipes derived), long double/float precision
  options removed, `-ascii` output option, `make all_no_gmt` without
  GMT libraries, references and citations in README and code.
- 2025: version 1.1.2.

## Copyright

Copyright (c) 2004-2026 Thorsten Becker, twb@ig.utexas.edu. All of
the HC code is distributed under the GNU General Public License, see
`COPYING`.

Parts of this software are derived from Numerical Recipes code: the
Gauss-Legendre quadrature (`rick_gauleg` in `src/rick_sh_c.c`) and the
FFT routines (`src/rick_fft_c.c`) used by the spherical harmonics
transforms, and, in the related tools, `src/fitxyee_util.c` and
`src/correl_nr_util.c`, which are copies of Numerical Recipes
routines. Copyright for those is with the Numerical Recipes authors.
(The 3 x 3 solver in `src/hc_matrix.c` used to be adapted from
Numerical Recipes as well, it is now an independent implementation.)
`src/ggrd_grdtrack_util.c` contains code derived from GMT, copyright
P. Wessel and W.H.F. Smith.
