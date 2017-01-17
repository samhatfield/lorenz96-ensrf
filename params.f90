module params
    use, intrinsic :: iso_fortran_env

    implicit none
    
    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    ! Frequency with which to write output
    integer, parameter :: write_freq = 50

    ! Output filename
    character(len=20) :: outfile = 'output.nc'

    ! Write values of X variables in output file?
    logical, parameter :: output_x = .true.

    ! Write values of Y variables in output file?
    logical, parameter :: output_y = .false.

    !Start and end times, and number of time steps
    real(dp), parameter :: dt = 0.005_dp
    real(dp), parameter :: fin = 12_dp
    integer, parameter :: n_steps = fin / dt

    ! State vector properties
    integer, parameter :: n_x = 8
    integer, parameter :: n_y = 8
    integer, parameter :: n_z = 8
    
    ! State vector dimension of full model
    integer, parameter :: truth_dim = n_x + n_x*n_y + n_x*n_y*n_z
    
    ! State vector dimension of parametrised model
    integer, parameter :: state_dim = n_x + n_x*n_y
    
    ! Number of ensemble members
    integer, parameter :: n_ens = 40

    ! Frequency of assimilations, i.e. 1 = every timestep
    integer, parameter :: assim_freq = 1

    ! Observation error variance
    real(dp), parameter :: y_var = 0.1_dp

    ! 1 = observe all Y variables, 2 = observe every other, 4 = observe every 4th
    integer, parameter :: y_skip = 1

    ! Dimension of observation vector
    integer, parameter :: obs_dim = n_x*n_y/y_skip

    ! Default significand
    integer, parameter :: sbits = 10

    ! Covariance inflation factor
    real(dp), parameter :: rho = 1.02_dp

    ! Covariance localisation length scales
    real(dp), parameter :: loc = 1.0_dp
    real(dp), parameter :: loc_y = loc*26
    real(dp), parameter :: loc_x = loc*8
end module params
