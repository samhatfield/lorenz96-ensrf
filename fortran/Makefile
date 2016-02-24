all:
	gfortran -c params.f90
	gfortran -c utils.f90
	gfortran -c lorenz96.f90
	gfortran -c metadata.f90
	gfortran -c *.f90 -Irpe/modules
	gfortran -o main *.o lapack/dgetrf/*.f lapack/dgetri/*.f -lblas -Lrpe/lib -lrpe 
	rm *.o
	rm *.mod

prof:
	gfortran -c -pg params.f90
	gfortran -c -pg utils.f90
	gfortran -c -pg lorenz96.f90
	gfortran -c -pg metadata.f90
	gfortran -c -pg *.f90 -Irpe/modules
	gfortran -o prof_main -pg *.o lapack/dgetrf/*.f lapack/dgetri/*.f -lblas -Lrpe/lib -lrpe 
	rm *.o
	rm *.mod

