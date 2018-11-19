# -*- Makefile -*-

FC      = ifort
OPTIM   = -g -C

CFITSIO = -L/mn/stornext/u3/hke/local/lib -lcfitsio
LAPACK  = -L/mn/stornext/u3/hke/local/lib -llapack -lblas
HEALPIX = -L/mn/stornext/u3/hke/local/lib -lhealpix
HEALINC = -I/mn/stornext/u3/hke/local/include
OUTPUT  = fittemps

OBJS    = fittemps.o

temp_fit: $(OBJS)
	$(FC) $(OBJS) $(HEALPIX) $(CFITSIO) -o $(OUTPUT)

# Compilation stage
%.o : %.f90
	$(FC) $(OPTIM) $(HEALINC) $(LAPACK) $(CFITSIO) -c $<

# Cleaning command
.PHONY: clean
clean:
	rm *.o *~ fittemps
