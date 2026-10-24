FC = mpif90
FFLAGS ?= -O2
ARITH ?= cmplx
LIB_DMRG ?= dmrg
DMRG_CFLAGS = $(shell pkg-config --cflags $(LIB_DMRG))
DMRG_LIBS = $(shell pkg-config --libs $(LIB_DMRG))
ifeq ($(ARITH),cmplx)
ARITH_FLAGS = -D_CMPLX
else ifeq ($(ARITH),dble)
ARITH_FLAGS =
else
$(error ARITH must be cmplx or dble, matching the selected nssDMRG library)
endif

.PHONY: all clean
all: atomic_bhz

atomic_bhz: atomic_bhz.f90 Makefile
	$(FC) $(FFLAGS) -cpp $(ARITH_FLAGS) $(DMRG_CFLAGS) $< $(DMRG_LIBS) -o $@

clean:
	rm -f atomic_bhz
