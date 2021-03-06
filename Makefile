# Default compiler
FC = gfortran

# Chooses full precision or reduced precision based on chosen target
double : PRECISION = real(dp)
single : PRECISION = real(sp)
rpe : PRECISION = type(rpe_var)

# Git revisions for L96-EnSRF and RPE emulator
GIT_REV_MAIN := $(shell git rev-parse HEAD)
GIT_REV_RPE := $(shell ( cd rpe ; git rev-parse HEAD ))

# Double precision - default
double: main

# Single precision
single: main

# Reduced precision
rpe: main

# Main target: main executable
main: main.o params.o lorenz96.o analysis.o setup.o utils.o observation.o io.o
	$(FC) $(COMPARGS) -o $@ $^ -lblas -lrpe -Lrpe/lib -L/usr/lib -lnetcdff -lnetcdf

# Dependencies
main.o: params.o lorenz96.o utils.o analysis.o setup.o observation.o io.o
lorenz96.o: params.o utils.o
utils.o: params.o
analysis.o: params.o utils.o observation.o
setup.o: params.o lorenz96.o utils.o
observation.o: params.o
io.o: params.o io.f90
	$(FC) $(COMPARGS) -c -cpp -DPRECISION='$(PRECISION)' -DGIT_REV_MAIN='"$(GIT_REV_MAIN)"'\
		-DPREC_STR='"$(PRECISION)"'	-DGIT_REV_RPE='"$(GIT_REV_RPE)"' io.f90 -o io.o\
		-Irpe/modules -I/usr/include


# Build rules
%.o: %.f90
	python parse.py $<
	$(FC) $(COMPARGS) -c -cpp -DPRECISION='$(PRECISION)' out.$< -o $(basename $<).o -Irpe/modules

.PHONY: clean
clean:
	rm -f *.o *.mod out.*.f90 main
