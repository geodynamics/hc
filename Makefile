#CFLAGS = -g
#
#
# makefile for experimental Hager & O'Connell routines
# and hcplates Ricard/Vigny type plate velocity inversions
#
# see source files for comments and reference to original authors
# 

#
# EDIT HERE FOR GMT VERSION 
#
include Makefile.include
#
#
#
#
ifdef ARCH
suffix=/$(ARCH)
endif
ifdef HC_HOME
prefix=$(HC_HOME)/
endif
# object file directory
ODIR = $(prefix)objects$(suffix)
#
# binary directory
BDIR = $(prefix)bin$(suffix)

# include files
OINCS = hc.h hc_filenames.h sh.h hc_constants.h
#
# Healpix stuff, comment out all if not wanted
#
# include flags
#HEAL_INC_DIR = $(HOME)/progs/src/Healpix_1.20/include/
#HEAL_INC_FLAGS = -I$(HEAL_INC_DIR)
#HEAL_LIBS = $(HOME)/progs/lib/$(ARCH)/libchealpix.a \
#		$(HOME)/progs/lib/$(ARCH)/libhealpix.a
#HEAL_LIB_FLAGS = -L/usr/local/src/cfitsio/lib/ -L/opt/cfitsio/lib/
#HEAL_LIBS_LINKLINE = -lchealpix -lhealpix -lcfitsio 
#HEAL_INCS = $(HEAL_INC_DIR)/myhealpix.h $(HEAL_INC_DIR)/chealpix.h 
#HEAL_DEFINES = -DHC_USE_HEALPIX
#
# Rick spherical harmonics stuff
#
#RICK_SRCS = rick_sh.f90 rick_fft.f90 rick_sh_c.c rick_fft_c.c
# new C version
RICK_SRCS = rick_sh_c.c rick_fft_c.c
#RICK_OBJS = $(ODIR)/rick_sh.o $(ODIR)/rick_sh_c.o  $(ODIR)/rick_fft.o $(ODIR)/rick_fft_c.o
#
#
#RICK_DEFINES = -DSH_RICK_DOUBLE_PRECISION 
# if -DNO_RICK_FORTRAN is defined, will only use C routines
RICK_DEFINES =  -DNO_RICK_FORTRAN 

RICK_OBJS = $(ODIR)/rick_sh_c.o $(ODIR)/rick_fft_c.o
RICK_OBJS_DBG = $(ODIR)/rick_sh_c.dbg.o $(ODIR)/rick_fft_c.dbg.o
RICK_INC_FLAGS = -I. 
RICK_INCS =  sh_rick_ftrn.h  sh_rick.h
RICK_LIB = $(ODIR)/librick.a $(ODIR)/librick.dbg.a

#
# PREM STUFF
#
PREM_SRCS = prem_util.c prem2r.c prem2rp.c prem2dsm.c 
PREM_OBJS = $(ODIR)/prem_util.o
PREM_DBG_OBJS = $(ODIR)/prem_util.dbg.o
# default PREM model file
PREM_DEFINES = -DPREM_MODEL_FILE=\"prem/prem.dat\"
PREM_INCS = prem.h
#
# GMT grd handling, now includes PREM stuff
#
GGRD_SRCS = ggrd_velinterpol.c ggrd_readgrds.c ggrd_grdtrack_util.c sh_exp_ggrd.c 
GGRD_OBJS = $(ODIR)/ggrd_velinterpol.o $(ODIR)/ggrd_readgrds.o \
	$(ODIR)/ggrd_grdtrack_util.o $(ODIR)/sh_exp_ggrd.o $(PREM_OBJS)

GGRD_OBJS_DBG = $(ODIR)/ggrd_velinterpol.dbg.o $(ODIR)/ggrd_readgrds.dbg.o $(ODIR)/ggrd_grdtrack_util.dbg.o 
GGRD_DEFINES = -I$(GMT4HOME)/include -I$(NETCDFHOME)/include  
GGRD_LIB_FLAGS = -L$(GMT4HOME)/lib -L$(NETCDFHOME)/lib 
GGRD_LIBS = $(ODIR)/libggrd.a $(ODIR)/libggrd.dfast.a $(ODIR)/libggrd.dbg.a
GGRD_LIBS_LINKLINE = $(GGRD_LIB_FLAGS) -lggrd -lpsl -lgmt -lnetcdf
GGRD_INCS = $(PREM_INCS)  ggrd_grdtrack_util.h ggrd_base.h ggrd_struc.h hc_ggrd_auto_proto.h

#
# other sources 
O_SOURCES = gaussp.c  print_gauss_lat.c rotvec2vel.c sh_corr.c simple_test.c  \
	spherepack_sh.c test_fft.c sh_test.c sh_power.c sh_model.c shana_sh.c
#
#
# Hager & O'Connell code
#
#
# C sources of subroutines (not main)
#
HC_SOURCES = sh_exp.c sh_model.c hc_init.c hc_solve.c hc_propagator.c \
	hc_polsol.c hc_matrix.c hc_torsol.c hc_output.c hc_input.c \
	hc_misc.c hc_extract_sh_layer.c  hc_extract_spatial.c hc_visc_scan.c

# all C sources without ggrd
C_SOURCES = $(HC_SOURCES) $(RICK_SRCS) $(O_SOURCES)
#
#
# objects for HC library
#
HC_OBJS = $(ODIR)/sh_exp.o $(ODIR)/sh_model.o $(ODIR)/hc_input.o \
	$(ODIR)/hc_polsol.o $(ODIR)/hc_matrix.o $(ODIR)/hc_torsol.o \
	$(ODIR)/hc_misc.o $(ODIR)/hc_init.o $(ODIR)/hc_propagator.o \
	$(ODIR)/hc_output.o $(ODIR)/hc_solve.o 

HC_OBJS_DBG = $(ODIR)/sh_exp.dbg.o $(ODIR)/sh_model.dbg.o $(ODIR)/hc_input.dbg.o \
	$(ODIR)/hc_polsol.dbg.o $(ODIR)/hc_matrix.dbg.o $(ODIR)/hc_torsol.dbg.o \
	$(ODIR)/hc_misc.dbg.o $(ODIR)/hc_init.dbg.o $(ODIR)/hc_propagator.dbg.o \
	$(ODIR)/hc_output.dbg.o $(ODIR)/hc_solve.dbg.o 

# HC libraries
HC_LIBS = $(ODIR)/libhc.a 
HC_LIBS_DEBUG =  $(ODIR)/libhc.dbg.a

LIB_FLAGS = $(HEAL_LIB_FLAGS) $(RICK_LIB_FLAGS) -L$(ODIR)/

#
INC_FLAGS =  $(HEAL_INC_FLAGS) $(ADD_FLAGS) $(RICK_INC_FLAGS)
#
# includes 
INCS = hc_auto_proto.h $(HEAL_INCS) $(RICK_INCS) $(OINCS)
#
# defines
DEFINES = $(RICK_DEFINES) $(HEAL_DEFINES) 
#
# libraries without ggrd
LIBS = $(HC_LIBS)  $(HEAL_LIBS) $(RICK_LIB)


all: all_no_gmt sh_tools hc_tools gmt_tools

all_no_gmt: $(ODIR) $(BDIR) libs hc_tools 

sh_tools: 	$(BDIR)/sh_syn $(BDIR)/sh_corr $(BDIR)/sh_ana $(BDIR)/sh_power \
	 $(BDIR)/sh_extract_layer $(BDIR)/gaussp

hc_tools: $(BDIR)/hc  $(BDIR)/hc_visc_scan $(BDIR)/hc_invert_dtopo \
	$(BDIR)/hc_extract_sh_layer  $(BDIR)/hc_extract_spatial \
	$(BDIR)/rotvec2vel $(BDIR)/print_gauss_lat

gmt_tools: $(BDIR)/grdgrd2correlation


weird_tools: $(BDIR)/convert_bernhard_dens

libs: $(ODIR) $(BDIR) hc_lib  $(HEAL_LIBS) $(RICK_LIB)

hc_lib: $(HC_LIBS)

ggrd_lib: ggrd_proto $(GGRD_LIBS)  

debug_libs: $(HC_LIBS_DEBUG)

really_all: proto all debug_libs weird_tools $(BDIR)/hc.dbg \
	hcplates $(BDIR)/ggrd_test $(BDIR)/grdinttester \
	$(BDIR)/prem2dsm $(BDIR)/prem2r  



proto: hc_auto_proto.h

ggrd_proto: hc_ggrd_auto_proto.h

hcplates: 
	cd hcplates; \
	make ;\
	cd ..


$(BDIR)/sh_test: $(LIBS) $(INCS) $(ODIR)/sh_test.o
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_test.o \
		-o $(BDIR)/sh_test -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)

$(BDIR)/sh_syn: $(LIBS) $(GGRD_LIBS) $(INCS) $(GGRD_INCS) $(ODIR)/sh_syn.o
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_syn.o \
		-o $(BDIR)/sh_syn -lhc -lrick $(HEAL_LIBS_LINKLINE) $(GGRD_LIBS_LINKLINE) -lm $(LDFLAGS)

$(BDIR)/sh_corr: $(LIBS) $(INCS) $(ODIR)/sh_corr.o $(PREM_OBJS) 
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_corr.o $(PREM_OBJS) \
		-o $(BDIR)/sh_corr -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)

$(BDIR)/sh_power: $(LIBS) $(INCS) $(PREM_OBJS)  $(ODIR)/sh_power.o
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_power.o \
		-o $(BDIR)/sh_power $(PREM_OBJS) -lhc -lrick $(HEAL_LIBS_LINKLINE)  -lm $(LDFLAGS)

$(BDIR)/sh_ana: $(LIBS) $(GGRD_LIBS) $(PREM_OBJS) $(INCS) $(GGRD_INCS) $(ODIR)/sh_ana.o
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_ana.o \
		-o $(BDIR)/sh_ana -lhc -lrick $(HEAL_LIBS_LINKLINE) \
		$(GGRD_LIBS_LINKLINE)  -lm $(LDFLAGS)

$(BDIR)/gaussp: $(LIBS) $(INCS) $(PREM_OBJS) $(ODIR)/gaussp.o
	$(CC) $(LIB_FLAGS) $(PREM_OBJS) $(ODIR)/gaussp.o \
		-o $(BDIR)/gaussp -lhc -lrick   -lm $(LDFLAGS)

$(BDIR)/sh_extract_layer: $(LIBS) $(INCS) $(PREM_OBJS) $(ODIR)/sh_extract_layer.o
	$(CC) $(LIB_FLAGS) $(ODIR)/sh_extract_layer.o \
		$(PREM_OBJS) -o $(BDIR)/sh_extract_layer \
	-lhc -lrick $(HEAL_LIBS_LINKLINE)  -lm $(LDFLAGS)

$(BDIR)/print_gauss_lat: print_gauss_lat.c $(PREM_OBJS) 
	$(CC) $(CFLAGS) $(PREM_OBJS) print_gauss_lat.c -o $(BDIR)/print_gauss_lat -lm $(INC_FLAGS) \
	$(LIB_FLAGS)   -lrick -lhc $(LDFLAGS) -lm

$(BDIR)/convert_bernhard_dens: convert_bernhard_dens.c
	$(CC) $(CFLAGS) convert_bernhard_dens.c $(PREM_OBJS) -o $(BDIR)/convert_bernhard_dens -lm $(INC_FLAGS) \
	$(LIB_FLAGS)   -lrick -lhc $(LDFLAGS)

$(BDIR)/rotvec2vel: rotvec2vel.c $(PREM_OBJS) 
	$(CC) $(CFLAGS) rotvec2vel.c $(PREM_OBJS)  -o $(BDIR)/rotvec2vel  -lm $(LDFLAGS)

$(BDIR)/prem2dsm: prem2dsm.c $(PREM_OBJS)
	$(CC) $(PREM_DEFINES) prem2dsm.c $(PREM_OBJS) -o $(BDIR)/prem2dsm -lm  $(LDFLAGS) 

$(BDIR)/prem2r: prem2r.c $(PREM_OBJS)
	$(CC)  $(PREM_DEFINES) prem2r.c $(PREM_OBJS) -o $(BDIR)/prem2r -lm $(LDFLAGS) 


$(BDIR)/hc: $(LIBS) $(INCS) $(ODIR)/hc.o $(PREM_OBJS)
	$(CC) $(LIB_FLAGS) $(ODIR)/hc.o -o $(BDIR)/hc \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) $(PREM_OBJS) -lm $(LDFLAGS) 

$(BDIR)/hc_visc_scan: $(LIBS) $(INCS) $(ODIR)/hc_visc_scan.o $(PREM_OBJS) 
	$(CC) $(LIB_FLAGS) $(ODIR)/hc_visc_scan.o -o $(BDIR)/hc_visc_scan \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) $(PREM_OBJS)  -lm $(LDFLAGS) 

$(BDIR)/hc_invert_dtopo: $(LIBS) $(INCS) $(ODIR)/hc_invert_dtopo.o $(PREM_OBJS) 
	$(CC) $(LIB_FLAGS) $(ODIR)/hc_invert_dtopo.o -o $(BDIR)/hc_invert_dtopo \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) $(PREM_OBJS)  \
		 -lm $(LDFLAGS) 

$(BDIR)/hc.dbg: $(LIBS) $(INCS) $(ODIR)/hc.dbg.o $(PREM_OBJS) 
	$(CC) $(LIB_FLAGS) $(ODIR)/hc.dbg.o -o $(BDIR)/hc.dbg \
		$(PREM_OBJS)	-lhc.dbg -lrick.dbg $(HEAL_LIBS_LINKLINE)   \
		 -lm $(LDFLAGS) 



$(BDIR)/test_fft: $(LIBS) $(INCS) $(ODIR)/test_fft.o
	$(CC) $(LIB_FLAGS) $(ODIR)/test_fft.o -o $(BDIR)/test_fft \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS) 

$(BDIR)/ggrd_test: $(LIBS) $(GGRD_LIBS) $(GGRD_INCS) $(INCS) ggrd_test.c
	$(CC) $(LIB_FLAGS) ggrd_test.c -o $(BDIR)/ggrd_test $(GGRD_DEFINES) \
		$(GGRD_LIBS_LINKLINE) -lhc -lrick -lm $(LDFLAGS) 

$(BDIR)/grdinttester: $(GGRD_INCS) $(LIBS) $(GGRD_LIBS) $(INCS) grdinttester.c
	$(CC) $(LIB_FLAGS) grdinttester.c $(GGRD_DEFINES) -o $(BDIR)/grdinttester \
		$(GGRD_LIBS_LINKLINE) -lhc -lrick -lm $(LDFLAGS) 

$(BDIR)/grdgrd2correlation: $(GGRD_INCS) $(LIBS) $(GGRD_LIBS) $(INCS) grdgrd2correlation.c fitxyee_util.c
	$(CC) $(LIB_FLAGS) grdgrd2correlation.c  fitxyee_util.c $(GGRD_DEFINES) -o $(BDIR)/grdgrd2correlation \
		$(GGRD_LIBS_LINKLINE) -lhc -lrick -lm $(LDFLAGS) 

$(BDIR)/hc_extract_sh_layer: $(LIBS) $(INCS) $(PREM_OBJS)  $(ODIR)/hc_extract_sh_layer.o
	$(CC) $(LIB_FLAGS) $(ODIR)/hc_extract_sh_layer.o  \
		-o $(BDIR)/hc_extract_sh_layer $(PREM_OBJS) \
		-lhc -lrick $(HEAL_LIBS_LINKLINE)  -lm $(LDFLAGS) 

$(BDIR)/hc_extract_spatial: $(LIBS) $(INCS) $(PREM_OBJS)  $(ODIR)/hc_extract_spatial.o
	$(CC) $(LIB_FLAGS) $(ODIR)/hc_extract_spatial.o $(PREM_OBJS)  \
		-o $(BDIR)/hc_extract_spatial \
		-lhc -lrick $(HEAL_LIBS_LINKLINE)   -lm $(LDFLAGS) 


#
# C function prototyper, strip out GMT version dependent things, 
# those are handled in other header
#
hc_auto_proto.h: 
	cproto  $(INC_FLAGS) $(DEFINES) -DGENERATE_PROTO  -f 2 -p $(C_SOURCES)  | \
		grep -v "void main("  | grep -v "int main(" > hc_auto_proto.h

hc_ggrd_auto_proto.h: 
	cproto  $(INC_FLAGS) $(GGRD_DEFINES) $(DEFINES) -DGENERATE_PROTO  -f 2 -p $(GGRD_SRCS)  | \
		grep -v "void main(" | grep -v "int main("   > hc_ggrd_auto_proto.h


$(ODIR):
	mkdir -p $(ODIR);

$(BDIR):
	mkdir -p $(BDIR);

clean:
	rm -f $(ODIR)/*.o  $(ODIR)/*.a $(BDIR)/*

super_clean:
	rm -f hc_auto_proto.h hc_ggrd_auto_proto.h $(ODIR)/*.o  $(ODIR)/*.a $(BDIR)/*	
#
# library
#

$(ODIR)/libhc.a: $(HC_OBJS) 
	$(AR) rv $(ODIR)/libhc.a $(HC_OBJS) 
	ranlib $(ODIR)/libhc.a


$(ODIR)/libhc.dbg.a: $(HC_OBJS_DBG)  $(PREM_DBG_OBJS)
	$(AR) rv $(ODIR)/libhc.dbg.a $(HC_OBJS_DBG)  $(PREM_DBG_OBJS)
	ranlib $(ODIR)/libhc.dbg.a


$(ODIR)/librick.a: $(RICK_OBJS)
	$(AR) rv $(ODIR)/librick.a $(RICK_OBJS)
	ranlib $(ODIR)/librick.a


$(ODIR)/librick.dbg.a: $(RICK_OBJS_DBG)
	$(AR) rv $(ODIR)/librick.dbg.a $(RICK_OBJS_DBG)
	ranlib $(ODIR)/librick.dbg.a

$(ODIR)/libggrd.a: $(GGRD_OBJS)
	$(AR) rv $(ODIR)/libggrd.a $(GGRD_OBJS)
	ranlib $(ODIR)/libggrd.a

$(ODIR)/libggrd.dfast.a: $(GGRD_OBJS)
	$(AR) rv $(ODIR)/libggrd.dfast.a $(GGRD_OBJS)
	ranlib $(ODIR)/libggrd.dfast.a


$(ODIR)/libggrd.dbg.a: $(GGRD_OBJS_DBG)
	$(AR) rv $(ODIR)/libggrd.dbg.a $(GGRD_OBJS_DBG)
	ranlib $(ODIR)/libggrd.dbg.a 

#
# object rules
#
$(ODIR)/hc_init.o: hc_init.c  $(INCS) $(PREM_INC)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(PREM_DEFINES) -c hc_init.c -o $(ODIR)/hc_init.o

$(ODIR)/hc_init.dbg.o: hc_init.c  $(INCS) $(PREM_INC)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) \
	$(DEFINES) $(PREM_DEFINES) -c hc_init.c -o $(ODIR)/hc_init.dbg.o

#
# ggrd/GMT stuff
#
$(ODIR)/ggrd_velinterpol.o: ggrd_velinterpol.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_velinterpol.c -o $(ODIR)/ggrd_velinterpol.o
$(ODIR)/ggrd_readgrds.o: ggrd_readgrds.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_readgrds.c -o $(ODIR)/ggrd_readgrds.o
$(ODIR)/ggrd_grdtrack_util.o: ggrd_grdtrack_util.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_grdtrack_util.c -o $(ODIR)/ggrd_grdtrack_util.o
$(ODIR)/sh_exp_ggrd.o: sh_exp_ggrd.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c sh_exp_ggrd.c -o $(ODIR)/sh_exp_ggrd.o
$(ODIR)/sh_ana.o: sh_ana.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c sh_ana.c -o $(ODIR)/sh_ana.o
$(ODIR)/sh_syn.o: sh_syn.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES)  \
	-c sh_syn.c -o $(ODIR)/sh_syn.o
# debug version
$(ODIR)/ggrd_velinterpol.dbg.o: ggrd_velinterpol.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_velinterpol.c -o $(ODIR)/ggrd_velinterpol.dbg.o
$(ODIR)/ggrd_readgrds.dbg.o: ggrd_readgrds.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_readgrds.c -o $(ODIR)/ggrd_readgrds.dbg.o
$(ODIR)/ggrd_grdtrack_util.dbg.o: ggrd_grdtrack_util.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c ggrd_grdtrack_util.c -o $(ODIR)/ggrd_grdtrack_util.dbg.o
$(ODIR)/sh_exp_ggrd.dbg.o: sh_exp_ggrd.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c sh_exp_ggrd.c -o $(ODIR)/sh_exp_ggrd.dbg.o
$(ODIR)/sh_ana.dbg.o: sh_ana.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES) \
	-c sh_ana.c -o $(ODIR)/sh_ana.dbg.o
$(ODIR)/sh_syn.dbg.o: sh_syn.c  $(INCS) $(GGRD_INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) $(GGRD_DEFINES)  \
	-c sh_syn.c -o $(ODIR)/sh_syn.dbg.o

#
# general rules
#
# debugging objects
$(ODIR)/%.dbg.o: %.c  $(INCS)
	$(CC) $(CFLAGS_DEBUG) -DHC_DEBUG $(INC_FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.f90 $(INCS)
	$(F90) $(F90FLAGS_DEBUG) -DHC_DEBUG $(DEFINES) -c $< -o $(ODIR)/$*.dbg.o


$(ODIR)/%.o: %.c  $(INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.f90 $(INCS)
	$(F90) $(F90FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.o

