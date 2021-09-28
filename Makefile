#TO BE CHANGED BY USER:
#$ DRIVER NAME without .f90 extension
#$ COMPILER: supported compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default if $HOME/.bin in the PATH)
EXE=ed_hm_bethe
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

define colorecho	
	@tput setaf 6
	@tput bold
	@echo $1
	@tput sgr0
endef



define lib_test
	_TEST+=$(shell pkg-config --cflags $(1);echo $$?)
	_TEST+=$(shell pkg-config --libs $(1);echo $$?)
	ifeq ($(_TEST),0)
	GLOB_INC+=$(shell pkg-config --cflags $(1))
	GLOB_LIB+=$(shell pkg-config --libs $(1))
	endif
endef



#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC:=$(shell pkg-config --cflags edlat dmft_ed slave_spins dmft_tools scifor)
GLOB_LIB:=$(shell pkg-config --libs edlat dmft_ed slave_spins dmft_tools scifor)
# $(call lib_test,"scifor")
# $(call lib_test,"dmft_tools")
# $(call lib_test,"slave_spins")
# $(call lib_test,"dmft_ed")
# # $(call lib_test,edlat)

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

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

all: FLAG:=${FFLAG} ${FPPMPI}
all:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

debug: FLAG:=${DFLAG} ${FPPMPI}
debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"



serial: FLAG:=${FFLAG} ${FPPSERIAL}
serial:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"

serial_debug: FLAG:=${DFLAG} ${FPPSERIAL}
serial_debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"





clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXE)



#########################################################################
