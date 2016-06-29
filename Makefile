# Default compiler
FC = gfortran

# Chooses full precision or reduced precision based on chosen target
double : PRECISION = real(dp)
single : PRECISION = real(sp)
rpe : PRECISION = type(rpe_var)

# Double precision - default
double: main

# Single precision
single: main

# Reduced precision
rpe: main

# Main target: main executable
main: main.o params.o lorenz96.o analysis.o metadata.o utils.o observation.o
	$(FC) -o $@ $^ -lblas -lrpe -Lrpe/lib

# Dependencies
main.o: params.o lorenz96.o utils.o analysis.o metadata.o observation.o
lorenz96.o: params.o utils.o
utils.o: params.o
analysis.o: params.o utils.o observation.o
metadata.o: params.o lorenz96.o
observation.o: params.o

# Build rules
%.o: %.f90
	python parse.py $<
	$(FC) -c -cpp -DPRECISION='$(PRECISION)' out.$< -o $(basename $<).o -Irpe/modules

.PHONY: clean
clean:
	rm -f *.o *.mod out.*.f90
