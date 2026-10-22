FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

# LIBRARIES TO BE INCLUDED
#$ LIB_TB: specify custom tight-binding library if any
#LIB_TB=honeytools
LIB_DMRG=dmrg

#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
ifdef LIB_DMRG
GLOB_INC+=$(shell pkg-config --cflags ${LIB_DMRG})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_DMRG})
endif

ifeq ($(PLAT),intel)
FFLAG=-O1 -ftz
DFLAG=-p -O0 -g
AFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif

ifeq ($(PLAT),gnu)
#FFLAG = -O1 -ffast-math -march=native -funroll-loops -ffree-line-length-none -g
FFLAG = -O3 -funroll-loops -ffree-line-length-none 
DFLAG = -w -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
AFLAG = -w -O0 -p -g  -fbacktrace -fwhole-file -fcheck=all -fbounds-check -fsanitize=address -fdebug-aux-vars -Wall -Waliasing -Wsurprising -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wdo-subscript -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wno-align-commons -Wno-overwrite-recursive -Wno-tabs -Wreal-q-constant -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wfrontend-loop-interchange -Wtarget-lifetime
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

define colorecho	
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef

all: bhz_post structure polarization


bhz_post: FLAG:=${FFLAG} ${FPPMPI}
bhz_post: energy corr chiPerp polarization structure wick
	$(call colorecho,"created drivers in $(DIREXE)", 2)
	@echo ""


energy: FLAG:=${FFLAG} ${FPPMPI}
energy:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/$@_bhz_1d $@_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""

corr: FLAG:=${FFLAG} ${FPPMPI}
corr:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/polarization_bhz_1d polarization_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""	

chiPerp: FLAG:=${FFLAG} ${FPPMPI}
chiPerp:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/polarization_bhz_1d polarization_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""

polarization: FLAG:=${FFLAG} ${FPPMPI}
polarization:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/$@_bhz_1d $@_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""

structure: FLAG:=${FFLAG} ${FPPMPI}
structure:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/$@_bhz_1d $@_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""

wick: FLAG:=${FFLAG} ${FPPMPI}
wick:
	$(call colorecho,"compiling $@ post-processing driver", 2)
	$(FC) -o $(DIREXE)/$@_bhz_1d $@_bhz_1d.f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo ""

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~



#########################################################################

# Standalone post-processing drivers (CMPLX nssDMRG).
.PHONY: structure polarization
