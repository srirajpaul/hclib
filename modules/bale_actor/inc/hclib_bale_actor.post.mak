# Post Makefile includes are the main part of a module's build system, allowing
# it to add flags to the overall project compile and link flags.
HCLIB_CFLAGS+=-I$(BALE_INSTALL)/include -I$(HCLIB_ROOT)/../modules/bale_actor/inc
HCLIB_LDFLAGS+=-L$(BALE_INSTALL)/lib -L$(HCLIB_ROOT)/../modules/bale_actor/lib
#HCLIB_LDLIBS+=-lconvey -lhclib_bale_actor
#HCLIB_LDLIBS+=-lconvey

