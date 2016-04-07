module params
    implicit none
    
    integer, parameter :: dp = kind(0.d0)

    !Start and end times, and number of time steps
    real(dp), parameter :: dt = 0.002_dp
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
    integer, parameter :: n_ens = 500

    ! Frequency of assimilations, i.e. 1 = every timestep
    integer, parameter :: assim_freq = 1

    ! Observation error variance
    real(dp), parameter :: y_var = 0.1_dp

    ! Inflation factor for tuning assimilation (see analysis.f90)
    real(dp), parameter :: rho = 1.0_dp

    ! Observe only Y variables
    integer, parameter :: obs_dim = n_x*n_y

    ! Default significand
    integer, parameter :: sbits = 15
end module params
