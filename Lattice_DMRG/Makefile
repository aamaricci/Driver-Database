#TO BE CHANGED BY USER:
#$ DRIVER NAME without .f90 extension
#$ COMPILER: supported compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default is $HOME/.bin)

EXE=dmrg_hm_1d
#EXE=dmrg_spin_1d
#EXE=dmrg_bhz_1d
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

ifdef LIB_TB
GLOB_INC+=$(shell pkg-config --cflags ${LIB_TB})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_TB})
endif

GLOB_INC+=$(shell pkg-config --cflags scifor)# libtcmalloc libprofiler)
GLOB_LIB+=$(shell pkg-config --libs   scifor)# libtcmalloc libprofiler)


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


serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)



debug: FLAG:=${DFLAG} ${FPPSERIAL} -D_DEBUG
debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/debug_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	$(call colorecho,"created debug_$(EXE) in  $(DIREXE)", 1)


mpi: FLAG:=${FFLAG} ${FPPMPI}
mpi:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC)   -o $(DIREXE)/mpi_$(EXE) $(EXE).f90 $(FLAG) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created mpi_$(EXE) in  $(DIREXE)", 1)

mpi_debug: FLAG:=${DFLAG} ${FPPMPI} -D_DEBUG
mpi_debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/mpi_debug_$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created mpi_debug_$(EXE) in  $(DIREXE)", 1)



clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/*$(EXE)



#########################################################################
