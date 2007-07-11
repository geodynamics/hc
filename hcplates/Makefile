CC=gcc
INCS=hcplates.h
INCS1=hc_findplate.h
INCS2=hc_ptrot.h
LIBS=-lm
LIST = hcplates.c hcplates_init.c zero_arrays.c plates_propagator.c get_loads.c read_unitrots_get_coeff.c integrate_torques.c crlb_vector_plm.c solve_rot.c svdcmp.c svbksb.c crlb_residual.c hcplates_write_output.c
LIST2 = findplate.c
LIST3 = ptrot.c hc_ptrot_mem.c hc_ptrot_init.c

all: hcplates hc_findplate hc_ptrot

hcplates: $(INCS)  $(LIST)

	${CC} $(INCS) -o hcplates $(LIST) $(CFLAGS) $(LIBS)

hc_findplate: $(INCS1) $(LIST2)

	gcc -g $(INCS1) -o hc_findplate $(LIST2) 

hc_ptrot: $(INCS2) $(LIST3)

	gcc -g $(INCS2) -o hc_ptrot $(LIST3) 



