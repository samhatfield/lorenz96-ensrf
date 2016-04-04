all:
	gfortran -c params.f90
	gfortran -c utils.f90 -Irpe/modules
	gfortran -c lorenz96.f90
	gfortran -c metadata.f90
	gfortran -c *.f90 -Irpe/modules
	gfortran -c rp_lapack/BLAS-3.6.0/*.f -Irpe/modules
	gfortran -c rp_lapack/rgetrf/*.f rp_lapack/rgetri/*.f -Irpe/modules
	gfortran -o main *.o -Lrpe/lib -lrpe
	rm *.o *.mod
