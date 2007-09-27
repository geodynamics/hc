#
#
# makefile for experimental Hager & O'Connell routines
# and hcplates Ricard/Vigny type plate velocity inversions
#
# see source files for comments and reference to original authors
# 
# 
#
#
LD = $(CC)
CFLAGS = $(CFLAGS_DEBUG) 

#
# EDIT HERE FOR GMT VERSION 
#
#
# for GMT3.4.5, use the next two lines
#GGRD_INC_FLAGS = -I$(GMTHOME)/include -I$(NETCDFHOME)/include 
#GGRD_LIBS_LINKLINE = -lggrd -lgmt -lnetcdf
# 
# for GMT version >= 4.1.2, uncomment the next two lines
GGRD_INC_FLAGS = -I$(GMTHOME)/include -I$(NETCDFHOME)/include -DUSE_GMT4
GGRD_LIBS_LINKLINE = -lggrd -lgmt -lpsl -lnetcdf 
#
#
#
#
# object file directory
ODIR = ../hc/objects/$(ARCH)/
#
#
# binary directory
BDIR = ../hc/bin/$(ARCH)/

# include files
OINCS = hc.h hc_filenames.h sh.h 
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
PREM_SRCS = prem_util.c
PREM_OBJS = $(ODIR)/prem_util.o
# default PREM model file
PREM_DEFINES = -DPREM_MODEL_FILE=\"$(PWD)/prem/prem.dat\"
PREM_INCS = prem.h
#
# GMT grd handling, now includes PREM stuff
#
GGRD_SRCS = ggrd_velinterpol.c ggrd_readgrds.c ggrd_grdtrack_util.c \
	$(PREM_SRCS)
GGRD_OBJS = $(ODIR)/ggrd_velinterpol.o $(ODIR)/ggrd_readgrds.o $(ODIR)/ggrd_grdtrack_util.o \
	$(PREM_OBJS)
GGRD_OBJS_DBG = $(ODIR)/ggrd_velinterpol.dbg.o $(ODIR)/ggrd_readgrds.dbg.o $(ODIR)/ggrd_grdtrack_util.dbg.o \
	$(PREM_OBJS)
GGRD_DEFINES = -I$(GMTHOME)/include -I$(NETCDFHOME)/include  \
	$(PREM_DEFINES)
GGRD_LIB_FLAGS = -L$(GMTHOME)/lib -L$(NETCDFHOME)/lib 
GGRD_LIBS = $(ODIR)/libggrd.a $(ODIR)/libggrd.dfast.a $(ODIR)/libggrd.dbg.a 
GGRD_INCS = $(PREM_INCS)  ggrd_grdtrack_util.h ggrd.h ggrd_struc.h

#
#
#
# Hager & O'Connell code
#
#
# C sources of subroutines (not main)
#
HC_SOURCES = sh_exp.c sh_model.c hc_init.c hc_solve.c hc_propagator.c \
	hc_polsol.c hc_matrix.c hc_torsol.c hc_output.c hc_input.c \
	hc_misc.c hc_extract_sh_layer.c 

# all C sources
C_SOURCES = $(HC_SOURCES) $(RICK_SRCS) $(GGRD_SRCS)
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

LIB_FLAGS = -L$(HOME)/progs/lib/$(ARCH)/ \
	$(HEAL_LIB_FLAGS) $(RICK_LIB_FLAGS) \
	$(GGRD_LIB_FLAGS) \
	-L$(ODIR)/

#
INC_FLAGS = -I$(HOME)/progs/include/  $(HEAL_INC_FLAGS) \
	$(RICK_INC_FLAGS) $(GGRD_INC_FLAGS) 
#
# includes 
INCS = hc_auto_proto.h $(HEAL_INCS) $(RICK_INCS)  $(GGRD_INCS) $(OINCS)
#
# defines
DEFINES = $(RICK_DEFINES) $(HEAL_DEFINES)  $(GGRD_DEFINES)
#
# libraries
LIBS = $(HC_LIBS) $(GGRD_LIBS) $(HEAL_LIBS) $(RICK_LIB)


all: dirs libs hc  hc_extract_sh_layer \
	sh_syn sh_ana sh_power 

libs: dirs hc_lib  $(HEAL_LIBS) $(RICK_LIB)

hc_lib: $(HC_LIBS) $(GGRD_LIBS)  

debug_libs: $(HC_LIBS_DEBUG)

really_all: proto all debug_libs hc.dbg hcplates ggrd_test grdinttester


proto: hc_auto_proto.h

hcplates: 
	cd hcplates; \
	make ;\
	cd ..


sh_test: $(LIBS) $(INCS) $(ODIR)/sh_test.o
	$(LD) $(LIB_FLAGS) $(ODIR)/sh_test.o \
		-o $(BDIR)/sh_test -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)

sh_syn: $(LIBS) $(INCS) $(ODIR)/sh_syn.o
	$(LD) $(LIB_FLAGS) $(ODIR)/sh_syn.o \
		-o $(BDIR)/sh_syn -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)

sh_power: $(LIBS) $(INCS) $(ODIR)/sh_power.o
	$(LD) $(LIB_FLAGS) $(ODIR)/sh_power.o \
		-o $(BDIR)/sh_power -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)

sh_ana: $(LIBS) $(INCS) $(ODIR)/sh_ana.o
	$(LD) $(LIB_FLAGS) $(ODIR)/sh_ana.o \
		-o $(BDIR)/sh_ana -lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS)


hc: $(LIBS) $(INCS) $(ODIR)/main.o 
	$(LD) $(LIB_FLAGS) $(ODIR)/main.o -o $(BDIR)/hc \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS) 

hc.dbg: $(LIBS) $(INCS) $(ODIR)/main.dbg.o 
	$(LD) $(LIB_FLAGS) $(ODIR)/main.dbg.o -o $(BDIR)/hc.dbg \
		-lhc.dbg -lrick.dbg \
	$(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS) 

test_fft: $(LIBS) $(INCS) $(ODIR)/test_fft.o
	$(LD) $(LIB_FLAGS) $(ODIR)/test_fft.o -o $(BDIR)/test_fft \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS) 

ggrd_test: $(LIBS) $(INCS) $(ODIR)/ggrd_test.o
	$(LD) $(LIB_FLAGS) $(ODIR)/ggrd_test.o -o $(BDIR)/ggrd_test \
		$(GGRD_LIBS_LINKLINE) -lhc -lrick -lm $(LDFLAGS) 

grdinttester: $(LIBS) $(INCS) $(ODIR)/grdinttester.o
	$(LD) $(LIB_FLAGS) $(ODIR)/grdinttester.o -o $(BDIR)/grdinttester \
		$(GGRD_LIBS_LINKLINE) -lhc -lrick -lm $(LDFLAGS) 

hc_extract_sh_layer: $(LIBS) $(INCS) $(ODIR)/hc_extract_sh_layer.o
	$(LD) $(LIB_FLAGS) $(ODIR)/hc_extract_sh_layer.o \
		-o $(BDIR)/hc_extract_sh_layer \
		-lhc -lrick $(HEAL_LIBS_LINKLINE) -lm $(LDFLAGS) 

#
# C function prototyper, strip out GMT version dependent things, 
# those are handled in other header
#
hc_auto_proto.h: 
	cproto  $(INC_FLAGS) $(DEFINES) -DGENERATE_PROTO  -f 2 -p *.c  | \
		grep -v "void main("  | \
		grep -v "ggrd_gt_bcr_init_loc(" | \
		grep -v "ggrd_grdtrack_interpolate(" | \
		grep -v "ggrd_grdtrack_init(" | \
	grep -v "int main(" > hc_auto_proto.h

dirs:
	if [ ! -s ../hc/ ]; then\
		mkdir ../hc/;\
	fi;
	if [ ! -s ../hc/objects/ ]; then\
		mkdir ../hc/objects;\
	fi;
	if [ ! -s ../hc/objects/$(ARCH)/ ]; then\
		mkdir ../hc/objects/$(ARCH);\
	fi;
	if [ ! -s ../hc/bin/ ];then\
		mkdir ../hc/bin;\
	fi;\
	if [ ! -s ../hc/bin/$(ARCH) ];then\
		mkdir ../hc/bin/$(ARCH);\
	fi

clean:
	rm $(ODIR)/*.o 
#
# library
#

$(ODIR)/libhc.a: $(HC_OBJS)
	$(AR) rv $(ODIR)/libhc.a $(HC_OBJS)

$(ODIR)/libhc.dbg.a: $(HC_OBJS_DBG)
	$(AR) rv $(ODIR)/libhc.dbg.a $(HC_OBJS_DBG)

$(ODIR)/librick.a: $(RICK_OBJS)
	$(AR) rv $(ODIR)/librick.a $(RICK_OBJS)

$(ODIR)/librick.dbg.a: $(RICK_OBJS_DBG)
	$(AR) rv $(ODIR)/librick.dbg.a $(RICK_OBJS_DBG)

$(ODIR)/libggrd.a: $(GGRD_OBJS)
	$(AR) rv $(ODIR)/libggrd.a $(GGRD_OBJS)

$(ODIR)/libggrd.dfast.a: $(GGRD_OBJS)
	$(AR) rv $(ODIR)/libggrd.dfast.a $(GGRD_OBJS)

$(ODIR)/libggrd.dbg.a: $(GGRD_OBJS_DBG)
	$(AR) rv $(ODIR)/libggrd.dbg.a $(GGRD_OBJS_DBG)

#
# object rules
#
$(ODIR)/%.o: %.c  $(INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.f90 $(INCS)
	$(F90) $(F90FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.o

# debugging objects
$(ODIR)/%.dbg.o: %.c  $(INCS)
	$(CC) $(CFLAGS_DEBUG) $(INC_FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.f90 $(INCS)
	$(F90) $(F90FLAGS_DEBUG) $(DEFINES) -c $< -o $(ODIR)/$*.dbg.o
