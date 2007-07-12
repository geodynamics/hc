
# object file directory
ODIR = ../../hc/objects/hcplates/$(ARCH)/
#
#
# binary directory
BDIR = ../../hc/bin/$(ARCH)/

INCS1=hcplates.h
INCS2=hc_findplate.h
INCS2=hc_ptrot.h

INCS = $(INCS1) $(INCS2) $(INCS3)

LIBS=-lm

LIST = $(ODIR)/hcplates.o $(ODIR)/hcplates_init.o $(ODIR)/zero_arrays.o \
	$(ODIR)/plates_propagator.o $(ODIR)/get_loads.o \
	$(ODIR)/read_unitrots_get_coeff.o $(ODIR)/integrate_torques.o \
	$(ODIR)/crlb_vector_plm.o $(ODIR)/solve_rot.o $(ODIR)/svdcmp.o \
	$(ODIR)/svbksb.o $(ODIR)/crlb_residual.o $(ODIR)/hcplates_write_output.o


LIST2 = $(ODIR)/findplate.o

LIST3 = $(ODIR)/ptrot.o $(ODIR)/hc_ptrot_mem.o $(ODIR)/hc_ptrot_init.o

all: dirs bins

dirs:
	if [ ! -s ../../hc/ ]; then\
		mkdir ../../hc/;\
	fi;\
	if [ ! -s ../../hc/objects/ ]; then\
		mkdir ../../hc/objects;\
	fi;\
	if [ ! -s ../../hc/objects/hcplates/ ]; then\
		mkdir ../../hc/objects/hcplates/;\
	fi;\
	if [ ! -s ../../hc/objects/hcplates/$(ARCH)/ ]; then\
		mkdir ../../hc/objects/hcplates/$(ARCH);\
	fi;\
	if [ ! -s ../../hc/bin/ ];then\
		mkdir ../../hc/bin;\
	fi;\
	if [ ! -s ../../hc/bin/$(ARCH) ];then\
		mkdir ../../hc/bin/$(ARCH);\
	fi



bins: $(BDIR)/hcplates $(BDIR)/hc_findplate $(BDIR)/hc_ptrot

$(BDIR)/hcplates: $(LIST)
	$(CC) -o $(BDIR)/hcplates $(LIST) $(LIBS)

$(BDIR)/hc_findplate: $(LIST2)
	$(CC)  -o $(BDIR)/hc_findplate $(LIST2)

$(BDIR)/hc_ptrot: $(INCS2) $(LIST3)
	$(CC)  -o $(BDIR)/hc_ptrot $(LIST3) 

$(ODIR)/%.o: %.c  $(INCS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(DEFINES) -c $< -o $(ODIR)/$*.o


