FC = gfortran
CC = cc
#LIBS = -lm -Wl,-framework -Wl,veclib 
LIBS = -llapack -lblas
PROG = SlaveSpin_DFT_seed.exe
FFLAGS = -O3

OBJ = SlaveSpin_DFT_seed.f broyd_grad.f broyd2_grad.f

all: $(PROG)

clobber:
	rm -f *.exe

clean:clobber
	rm -f SlaveSpin_DFT_seed.exe

$(PROG): $(OBJ)
	$(FC) $(FFLAGS) -o $(PROG) $(OBJ) $(LIBS)


.f.o : 
	$(FC) $(FFLAGS) -c  $<

.c.o : 
	$(CC) -O2  -c  $<
