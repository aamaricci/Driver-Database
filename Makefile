#TO BE CHANGED BY USER:
#$ EXE: the DRIVER NAME 
#$ COMPILER: supported compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default if $HOME/.bin in the PATH)
EXE=solveBHZ_reticolo.f90
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

define colorecho	
	@tput setaf 6
	@tput bold
	@echo $1
	@tput sgr0
endef





#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC:=$(shell pkg-config --cflags dmft_ed slave_spins dmft_tools scifor)
GLOB_LIB:=$(shell pkg-config --libs dmft_ed slave_spins dmft_tools scifor)

ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif
ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O2 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI	
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

EXEC=$(shell basename -s .f90 ${EXE})

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

all: FLAG:=${FFLAG} ${FPPMPI}
all:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

debug: FLAG:=${DFLAG} ${FPPMPI}
debug:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"



serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

serial_debug: FLAG:=${DFLAG} ${FPPSERIAL}
serial_debug:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"





clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXEC)



#########################################################################
