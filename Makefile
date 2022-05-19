#TO BE CHANGED BY USER:
#$ EXE: the DRIVER NAME
#$ FC: the fortran compiler, default is mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default if $HOME/.bin)
EXE=1d_Kondo_chain.f90
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

# LIBRARIES TO BE INCLUDED
#$ LIB_ED: either use *edlat* or *dmft_ed* (until we fix the naming conventions)
#$ LIB_SS: specify slave spins library is any
LIB_ED=dmft_ed
LIB_SS=slave_spins





#########################################################################
#########################################################################
#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC=
GLOB_LIB=

ifdef LIB_ED
GLOB_INC+=$(shell pkg-config --cflags ${LIB_ED})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_ED})
endif

ifdef LIB_SS
GLOB_INC+=$(shell pkg-config --cflags ${LIB_SS})
GLOB_LIB+=$(shell pkg-config --libs ${LIB_SS})
endif

GLOB_INC+=$(shell pkg-config --cflags dmft_tools scifor)
GLOB_LIB+=$(shell pkg-config --libs   dmft_tools scifor)




ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPSERIAL =-fpp -D_
FPPMPI =-fpp -D_	
endif
ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O2 -p -g -fimplicit-none -Wsurprising -Wall -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI	
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

EXEC=$(shell basename ${EXE} .f90)

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

define colorecho	
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef



all: FLAG:=${OFLAG} ${FPPMPI}
all:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ", 2)
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXEC) in  $(DIREXE)", 1)

debug: FLAG:=${DFLAG} ${FPPMPI}
debug:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ", 2)
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXEC) in  $(DIREXE)", 1)



serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ", 2)
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXEC) in  $(DIREXE)", 1)

serial_debug: FLAG:=${DFLAG} ${FPPSERIAL}
serial_debug:
	@echo ""
	$(call colorecho,"compiling $(EXEC).f90 ", 2)
	@echo ""
	$(FC) $(FLAG) $(EXEC).f90 -o $(DIREXE)/$(EXEC) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXEC) in  $(DIREXE)", 1)




clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXEC)



#########################################################################
