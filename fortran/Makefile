all:
	gfortran -c params.f90
	gfortran -c utils.f90
	gfortran -c lorenz96.f90
	gfortran -c metadata.f90
	gfortran -c *.f90 -Irpe/modules
	gfortran -o main *.o lapack/dgetrf/*.f lapack/dgetri/*.f -lblas -Lrpe/lib -lrpe 
	rm *.o
	rm *.mod
