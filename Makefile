# Default compiler
FC = gfortran

# Chooses full precision or reduced precision based on chosen target
full : DOUBLE_OR_RPE = real(dp)
rpe : DOUBLE_OR_RPE = type(rpe_var)

# LAPACK files required for full precision matrix inverse
LAPACK = $(patsubst %.f,%.o,$(wildcard lapack/dgetrf/*.f lapack/dgetri/*.f))

# Reduced precision LAPACK and BLAS files required for reduced precision matrix inverse
RP_LAPACK = $(patsubst %.f,%.o,$(wildcard rp_lapack/rgetrf/*.f rp_lapack/rgetri/*.f))
BLAS = $(patsubst %.f,%.o,$(wildcard rp_lapack/BLAS-3.6.0/*.f))

# Full precision - default
full: main

# Reduced precision
rpe: main

# Main target: main executable
main: main.o params.o lorenz96.o analysis.o metadata.o utils.o observation.o\
	$(LAPACK) $(RP_LAPACK) $(BLAS)
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
	$(FC) -c -cpp -DDOUBLE_OR_RPE='$(DOUBLE_OR_RPE)' $< -Irpe/modules

rp_lapack/rgetr%.o: rp_lapack/rgetr%.f
	$(FC) -c -o $@ $< -Irpe/modules

rp_lapack/BLAS-3.6.0/%.o: rp_lapack/BLAS-3.6.0/%.f
	$(FC) -c -o $@ $< -Irpe/modules

.PHONY: clean
clean:
	rm -f *.o *.mod $(LAPACK) $(RP_LAPACK) $(BLAS)
