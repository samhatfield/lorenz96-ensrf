all:
	gfortran -c -Irpe/modules main.f90
	gfortran -o main main.o lapack/dgetrf/*.f lapack/dgetri/*.f -lblas -Lrpe/lib -lrpe 
